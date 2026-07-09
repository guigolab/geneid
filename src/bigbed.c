/*************************************************************************
*   Module: bigbed -- minimal random-access reader for UCSC bigBed (BBI).*
*                                                                        *
*   Only what geneid needs: open + header + chrom B+-tree, then a range  *
*   query that walks the R-tree index and decodes ONLY the data blocks   *
*   overlapping [start,end). Little-endian files only (what every writer *
*   emits); zlib-compressed blocks (uncompressBufSize>0) are inflated.   *
*                                                                        *
*   Format reference: Kent et al., "BigWig and BigBed", Bioinformatics   *
*   2010; the BBI on-disk layout (bbiFile.h / cirTree.h in kentUtils).   *
*   This file is part of the geneid Distribution.                        *
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <zlib.h>

#include "bigbed.h"

#define BIGBED_MAGIC 0x8789F2EBu
#define BPT_MAGIC    0x78CA8C91u   /* chrom B+-tree */
#define CIRTREE_MAGIC 0x2468ACE0u  /* R-tree index  */

typedef struct { char* name; unsigned id; unsigned size; } bbChrom;

struct BigBed {
  int   fd;
  unsigned uncompressBufSize;      /* 0 => data blocks are stored uncompressed */
  unsigned long fullIndexOffset;   /* start of the R-tree */
  bbChrom* chroms;
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
/* Recursively load every (name,id,size) leaf into bb->chroms. */
static int bptLoad(BigBed* bb, unsigned long off, unsigned keySize){
  unsigned char hdr[4];
  if (readAt(bb->fd, hdr, 4, off)) return -1;
  int isLeaf = hdr[0];
  int count  = rdU16(hdr + 2);
  unsigned long p = off + 4;
  if (isLeaf){
    unsigned char* item = malloc(keySize + 8);
    if (!item) return -1;
    for (int i = 0; i < count; i++){
      if (readAt(bb->fd, item, keySize + 8, p)){ free(item); return -1; }
      p += keySize + 8;
      bbChrom* c = &bb->chroms[bb->nChroms++];
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
      if (readAt(bb->fd, item, keySize + 8, p)){ free(item); free(kids); return -1; }
      p += keySize + 8;
      kids[i] = rdU64(item + keySize);   /* child node offset */
    }
    free(item);
    for (int i = 0; i < count; i++)
      if (bptLoad(bb, kids[i], keySize)){ free(kids); return -1; }
    free(kids);
  }
  return 0;
}

BigBed* bbOpen(const char* path){
  BigBed* bb = calloc(1, sizeof(BigBed));
  if (!bb) return NULL;
  bb->fd = open(path, O_RDONLY);
  if (bb->fd < 0){ free(bb); return NULL; }

  unsigned char hdr[64];
  if (readAt(bb->fd, hdr, 64, 0) || rdU32(hdr) != BIGBED_MAGIC){ bbClose(bb); return NULL; }
  unsigned long chromTreeOffset = rdU64(hdr + 8);
  bb->fullIndexOffset           = rdU64(hdr + 24);
  bb->uncompressBufSize         = rdU32(hdr + 52);

  /* chrom B+-tree header: magic(4) blockSize(4) keySize(4) valSize(4)
     itemCount(8) reserved(8) -> nodes follow at +32 */
  unsigned char bpt[32];
  if (readAt(bb->fd, bpt, 32, chromTreeOffset) || rdU32(bpt) != BPT_MAGIC){ bbClose(bb); return NULL; }
  unsigned keySize   = rdU32(bpt + 8);
  unsigned long nItems = rdU64(bpt + 16);
  bb->chroms = calloc(nItems ? nItems : 1, sizeof(bbChrom));
  if (!bb->chroms){ bbClose(bb); return NULL; }
  if (bptLoad(bb, chromTreeOffset + 32, keySize)){ bbClose(bb); return NULL; }
  return bb;
}

void bbClose(BigBed* bb){
  if (!bb) return;
  if (bb->fd >= 0) close(bb->fd);
  for (int i = 0; i < bb->nChroms; i++) free(bb->chroms[i].name);
  free(bb->chroms);
  free(bb);
}

static int chromId(BigBed* bb, const char* name, unsigned* id){
  for (int i = 0; i < bb->nChroms; i++)
    if (!strcmp(bb->chroms[i].name, name)){ *id = bb->chroms[i].id; return 1; }
  return 0;
}

/* Does block range [(sC,sB),(eC,eB)] overlap query (qC,[qs,qe])? Closed-interval
   overlap (a feature touching the query edge counts), matching the reference
   readers (pybigtools/kent) so boundary features at a split edge are not lost. */
static int cirOverlap(unsigned sC, unsigned sB, unsigned eC, unsigned eB,
                      unsigned qC, unsigned qs, unsigned qe){
  /* block entirely before the query? (ends strictly before qs) */
  if (eC < qC || (eC == qC && eB < qs)) return 0;
  /* block entirely after the query? (starts strictly after qe) */
  if (sC > qC || (sC == qC && sB > qe)) return 0;
  return 1;
}

/* Decode one data block (already read into `data`, length n) and report the
   records overlapping (qId,[qs,qe)). Returns #reported. */
static long emitBlock(const unsigned char* data, size_t n, unsigned qId,
                      unsigned qs, unsigned qe, bbRecordCB cb, void* ud){
  long count = 0;
  size_t i = 0;
  while (i + 12 <= n){
    unsigned cId = rdU32(data + i);
    unsigned s   = rdU32(data + i + 4);
    unsigned e   = rdU32(data + i + 8);
    i += 12;
    const char* rest = (const char*)(data + i);       /* NUL-terminated */
    size_t rlen = strnlen(rest, n - i);
    if (cId == qId && s <= qe && e >= qs) { cb((long)s, (long)e, rest, ud); count++; }
    i += rlen + 1;
  }
  return count;
}

/* Walk the R-tree from node `off`, decoding overlapping leaf blocks. */
static long cirQuery(BigBed* bb, unsigned long off, unsigned qId, unsigned qs,
                     unsigned qe, bbRecordCB cb, void* ud){
  unsigned char hdr[4];
  if (readAt(bb->fd, hdr, 4, off)) return -1;
  int isLeaf = hdr[0];
  int count  = rdU16(hdr + 2);
  unsigned long p = off + 4;
  long total = 0;

  if (isLeaf){
    unsigned char it[32];                              /* 4*4 + 8 + 8 */
    for (int i = 0; i < count; i++){
      if (readAt(bb->fd, it, 32, p)) return -1;
      p += 32;
      unsigned sC = rdU32(it), sB = rdU32(it+4), eC = rdU32(it+8), eB = rdU32(it+12);
      unsigned long dOff = rdU64(it+16), dSize = rdU64(it+24);
      if (!cirOverlap(sC, sB, eC, eB, qId, qs, qe)) continue;
      unsigned char* raw = malloc(dSize);
      if (!raw) return -1;
      if (readAt(bb->fd, raw, dSize, dOff)){ free(raw); return -1; }
      long r;
      if (bb->uncompressBufSize){
        uLongf out = bb->uncompressBufSize;
        unsigned char* ub = malloc(out);
        if (!ub){ free(raw); return -1; }
        if (uncompress(ub, &out, raw, dSize) != Z_OK){ free(ub); free(raw); return -1; }
        r = emitBlock(ub, out, qId, qs, qe, cb, ud);
        free(ub);
      } else {
        r = emitBlock(raw, dSize, qId, qs, qe, cb, ud);
      }
      free(raw);
      if (r < 0) return -1;
      total += r;
    }
  } else {
    unsigned char it[24];                              /* 4*4 + 8 */
    unsigned long* kids = malloc(sizeof(unsigned long) * count);
    unsigned* keep = malloc(sizeof(unsigned) * count);
    int nk = 0;
    for (int i = 0; i < count; i++){
      if (readAt(bb->fd, it, 24, p)){ free(kids); free(keep); return -1; }
      p += 24;
      unsigned sC = rdU32(it), sB = rdU32(it+4), eC = rdU32(it+8), eB = rdU32(it+12);
      if (cirOverlap(sC, sB, eC, eB, qId, qs, qe)){ kids[nk] = rdU64(it+16); keep[nk] = 1; nk++; }
    }
    (void)keep;
    for (int i = 0; i < nk; i++){
      long r = cirQuery(bb, kids[i], qId, qs, qe, cb, ud);
      if (r < 0){ free(kids); free(keep); return -1; }
      total += r;
    }
    free(kids); free(keep);
  }
  return total;
}

long bbQuery(BigBed* bb, const char* chrom, long start, long end,
             bbRecordCB cb, void* ud){
  unsigned id;
  if (!chromId(bb, chrom, &id)) return 0;
  /* R-tree header is 48 bytes; root node follows. */
  return cirQuery(bb, bb->fullIndexOffset + 48, id, (unsigned)start, (unsigned)end, cb, ud);
}
