/*************************************************************************
*   Module: bigwig -- minimal random-access reader for UCSC bigWig (BBI). *
*                                                                        *
*   Only what geneid needs: open + header + chrom B+-tree, then a range  *
*   query that walks the R-tree index and decodes ONLY the wiggle        *
*   sections overlapping [start,end). Little-endian files only (what      *
*   every writer emits); zlib-compressed blocks (uncompressBufSize>0)     *
*   are inflated. Handles the three section encodings a bigWig can hold:  *
*   bedGraph (1), varStep (2), fixedStep (3).                            *
*                                                                        *
*   The BBI container (64-byte header, chrom B+-tree, R-tree index) is    *
*   identical to bigBed; this reader mirrors src/bigbed.c and diverges    *
*   only in the leaf-block decoder (emitSection vs emitBlock). Kept as a  *
*   separate self-contained file to match bigbed.c; factoring the shared  *
*   container out is a possible future cleanup.                          *
*                                                                        *
*   Format reference: Kent et al., "BigWig and BigBed", Bioinformatics    *
*   2010; the BBI on-disk layout (bbiFile.h / cirTree.h in kentUtils).    *
*   This file is part of the geneid Distribution.                        *
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <zlib.h>

#include "bigwig.h"

#define BIGWIG_MAGIC  0x888FFC26u
#define BPT_MAGIC     0x78CA8C91u   /* chrom B+-tree */
#define CIRTREE_MAGIC 0x2468ACE0u   /* R-tree index  */

/* wiggle section types (byte at section-header offset 20) */
#define SECT_BEDGRAPH 1
#define SECT_VARSTEP  2
#define SECT_FIXEDSTEP 3

typedef struct { char* name; unsigned id; unsigned size; } bwChrom;

struct BigWig {
  int   fd;
  unsigned uncompressBufSize;      /* 0 => data blocks are stored uncompressed */
  unsigned long fullIndexOffset;   /* start of the R-tree */
  bwChrom* chroms;
  int   nChroms;
};

/* ---- little-endian scalar reads from a byte buffer -------------------------*/
static unsigned rdU16(const unsigned char* p){ return p[0] | (p[1]<<8); }
static unsigned rdU32(const unsigned char* p){
  return (unsigned)p[0] | ((unsigned)p[1]<<8) | ((unsigned)p[2]<<16) | ((unsigned)p[3]<<24);
}
static unsigned long rdU64(const unsigned char* p){
  return (unsigned long)rdU32(p) | ((unsigned long)rdU32(p+4) << 32);
}
static float rdF32(const unsigned char* p){        /* IEEE-754 little-endian */
  unsigned u = rdU32(p);
  float f;
  memcpy(&f, &u, 4);
  return f;
}

/* pread exactly n bytes at off; return 0 on success, -1 on short/failed read */
static int readAt(int fd, void* buf, size_t n, unsigned long off){
  size_t got = 0;
  while (got < n){
    ssize_t r = pread(fd, (char*)buf + got, n - got, (off_t)(off + got));
    if (r <= 0) return -1;
    got += (size_t)r;
  }
  return 0;
}

/* ---- chrom B+-tree -------------------------------------------------------- */
/* Recursively load every (name,id,size) leaf into bw->chroms. */
static int bptLoad(BigWig* bw, unsigned long off, unsigned keySize){
  unsigned char hdr[4];
  if (readAt(bw->fd, hdr, 4, off)) return -1;
  int isLeaf = hdr[0];
  int count  = rdU16(hdr + 2);
  unsigned long p = off + 4;
  if (isLeaf){
    unsigned char* item = malloc(keySize + 8);
    if (!item) return -1;
    for (int i = 0; i < count; i++){
      if (readAt(bw->fd, item, keySize + 8, p)){ free(item); return -1; }
      p += keySize + 8;
      bwChrom* c = &bw->chroms[bw->nChroms++];
      c->name = malloc(keySize + 1);
      memcpy(c->name, item, keySize);
      c->name[keySize] = '\0';           /* names are NUL-padded to keySize */
      c->id   = rdU32(item + keySize);
      c->size = rdU32(item + keySize + 4);
    }
    free(item);
  } else {
    unsigned char* item = malloc(keySize + 8);
    if (!item) return -1;
    unsigned long* kids = malloc(sizeof(unsigned long) * count);
    if (!kids){ free(item); return -1; }
    for (int i = 0; i < count; i++){
      if (readAt(bw->fd, item, keySize + 8, p)){ free(item); free(kids); return -1; }
      p += keySize + 8;
      kids[i] = rdU64(item + keySize);   /* child node offset */
    }
    free(item);
    for (int i = 0; i < count; i++)
      if (bptLoad(bw, kids[i], keySize)){ free(kids); return -1; }
    free(kids);
  }
  return 0;
}

BigWig* bwOpen(const char* path){
  BigWig* bw = calloc(1, sizeof(BigWig));
  if (!bw) return NULL;
  bw->fd = open(path, O_RDONLY);
  if (bw->fd < 0){ free(bw); return NULL; }

  unsigned char hdr[64];
  if (readAt(bw->fd, hdr, 64, 0) || rdU32(hdr) != BIGWIG_MAGIC){ bwClose(bw); return NULL; }
  unsigned long chromTreeOffset = rdU64(hdr + 8);
  bw->fullIndexOffset           = rdU64(hdr + 24);
  bw->uncompressBufSize         = rdU32(hdr + 52);

  /* chrom B+-tree header: magic(4) blockSize(4) keySize(4) valSize(4)
     itemCount(8) reserved(8) -> nodes follow at +32 */
  unsigned char bpt[32];
  if (readAt(bw->fd, bpt, 32, chromTreeOffset) || rdU32(bpt) != BPT_MAGIC){ bwClose(bw); return NULL; }
  unsigned keySize   = rdU32(bpt + 8);
  unsigned long nItems = rdU64(bpt + 16);
  bw->chroms = calloc(nItems ? nItems : 1, sizeof(bwChrom));
  if (!bw->chroms){ bwClose(bw); return NULL; }
  if (bptLoad(bw, chromTreeOffset + 32, keySize)){ bwClose(bw); return NULL; }
  return bw;
}

void bwClose(BigWig* bw){
  if (!bw) return;
  if (bw->fd >= 0) close(bw->fd);
  for (int i = 0; i < bw->nChroms; i++) free(bw->chroms[i].name);
  free(bw->chroms);
  free(bw);
}

static int chromId(BigWig* bw, const char* name, unsigned* id){
  for (int i = 0; i < bw->nChroms; i++)
    if (!strcmp(bw->chroms[i].name, name)){ *id = bw->chroms[i].id; return 1; }
  return 0;
}

/* Does block bound [(sC,sB),(eC,eB)] overlap query (qC,[qs,qe))? Conservative
   closed-interval test used only to PRUNE the R-tree walk -- never prunes a
   block that could hold a needed item; the exact half-open filter is applied
   per item in emitSection, so a loose prune here stays correct. */
static int cirOverlap(unsigned sC, unsigned sB, unsigned eC, unsigned eB,
                      unsigned qC, unsigned qs, unsigned qe){
  if (eC < qC || (eC == qC && eB < qs)) return 0;   /* entirely before query */
  if (sC > qC || (sC == qC && sB > qe)) return 0;   /* entirely after  query */
  return 1;
}

/* Report one interval if it overlaps (qId,[qs,qe)) half-open. */
static long emitInterval(unsigned cId, long s, long e, float v, unsigned qId,
                         unsigned qs, unsigned qe, bwIntervalCB cb, void* ud){
  if (cId == qId && s < (long)qe && e > (long)qs){ cb(s, e, v, ud); return 1; }
  return 0;
}

/* Decode one data block (already read into `data`, length n) -- a run of one or
   more wiggle sections -- and report the intervals overlapping (qId,[qs,qe)).
   Returns #reported, or -1 on a malformed section. */
static long emitSection(const unsigned char* data, size_t n, unsigned qId,
                        unsigned qs, unsigned qe, bwIntervalCB cb, void* ud){
  long count = 0;
  size_t i = 0;
  while (i + 24 <= n){                                /* section header = 24B */
    unsigned cId    = rdU32(data + i);
    unsigned cStart = rdU32(data + i + 4);
    /* data + i + 8  = chromEnd (unused: item coords are self-describing) */
    unsigned step   = rdU32(data + i + 12);
    unsigned span   = rdU32(data + i + 16);
    unsigned type   = data[i + 20];
    unsigned items  = rdU16(data + i + 22);
    i += 24;

    for (unsigned k = 0; k < items; k++){
      long s, e;
      float v;
      if (type == SECT_BEDGRAPH){
        if (i + 12 > n) return -1;
        s = (long)rdU32(data + i);
        e = (long)rdU32(data + i + 4);
        v = rdF32(data + i + 8);
        i += 12;
      } else if (type == SECT_VARSTEP){
        if (i + 8 > n) return -1;
        s = (long)rdU32(data + i);
        e = s + (long)span;
        v = rdF32(data + i + 4);
        i += 8;
      } else if (type == SECT_FIXEDSTEP){
        if (i + 4 > n) return -1;
        s = (long)cStart + (long)k * (long)step;
        e = s + (long)span;
        v = rdF32(data + i);
        i += 4;
      } else {
        return -1;                                    /* unknown section type */
      }
      count += emitInterval(cId, s, e, v, qId, qs, qe, cb, ud);
    }
  }
  return count;
}

/* Walk the R-tree from node `off`, decoding overlapping leaf blocks. */
static long cirQuery(BigWig* bw, unsigned long off, unsigned qId, unsigned qs,
                     unsigned qe, bwIntervalCB cb, void* ud){
  unsigned char hdr[4];
  if (readAt(bw->fd, hdr, 4, off)) return -1;
  int isLeaf = hdr[0];
  int count  = rdU16(hdr + 2);
  unsigned long p = off + 4;
  long total = 0;

  if (isLeaf){
    unsigned char it[32];                              /* 4*4 + 8 + 8 */
    for (int i = 0; i < count; i++){
      if (readAt(bw->fd, it, 32, p)) return -1;
      p += 32;
      unsigned sC = rdU32(it), sB = rdU32(it+4), eC = rdU32(it+8), eB = rdU32(it+12);
      unsigned long dOff = rdU64(it+16), dSize = rdU64(it+24);
      if (!cirOverlap(sC, sB, eC, eB, qId, qs, qe)) continue;
      unsigned char* raw = malloc(dSize);
      if (!raw) return -1;
      if (readAt(bw->fd, raw, dSize, dOff)){ free(raw); return -1; }
      long r;
      if (bw->uncompressBufSize){
        uLongf out = bw->uncompressBufSize;
        unsigned char* ub = malloc(out);
        if (!ub){ free(raw); return -1; }
        if (uncompress(ub, &out, raw, dSize) != Z_OK){ free(ub); free(raw); return -1; }
        r = emitSection(ub, out, qId, qs, qe, cb, ud);
        free(ub);
      } else {
        r = emitSection(raw, dSize, qId, qs, qe, cb, ud);
      }
      free(raw);
      if (r < 0) return -1;
      total += r;
    }
  } else {
    unsigned char it[24];                              /* 4*4 + 8 */
    unsigned long* kids = malloc(sizeof(unsigned long) * count);
    if (!kids) return -1;
    int nk = 0;
    for (int i = 0; i < count; i++){
      if (readAt(bw->fd, it, 24, p)){ free(kids); return -1; }
      p += 24;
      unsigned sC = rdU32(it), sB = rdU32(it+4), eC = rdU32(it+8), eB = rdU32(it+12);
      if (cirOverlap(sC, sB, eC, eB, qId, qs, qe)) kids[nk++] = rdU64(it+16);
    }
    for (int i = 0; i < nk; i++){
      long r = cirQuery(bw, kids[i], qId, qs, qe, cb, ud);
      if (r < 0){ free(kids); return -1; }
      total += r;
    }
    free(kids);
  }
  return total;
}

long bwQuery(BigWig* bw, const char* chrom, long start, long end,
             bwIntervalCB cb, void* ud){
  unsigned id;
  if (!chromId(bw, chrom, &id)) return 0;
  /* R-tree header is 48 bytes; root node follows. */
  return cirQuery(bw, bw->fullIndexOffset + 48, id, (unsigned)start, (unsigned)end, cb, ud);
}
