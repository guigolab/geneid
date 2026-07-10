/* Unit test for emitSection() -- the wiggle-section decoder inside bigwig.c.
   pybigtools only writes bedGraph sections, so the varStep/fixedStep branches
   (and their coordinate arithmetic) are not covered by the oracle test; this
   drives all three section encodings directly with hand-built byte buffers.

   Includes the translation unit so the static emitSection is reachable. */
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include "../../src/bigwig.c"

/* little-endian encoders into a growing buffer */
static void pU32(unsigned char* b, size_t* i, unsigned v){
  b[(*i)++]=v&0xff; b[(*i)++]=(v>>8)&0xff; b[(*i)++]=(v>>16)&0xff; b[(*i)++]=(v>>24)&0xff;
}
static void pU16(unsigned char* b, size_t* i, unsigned v){ b[(*i)++]=v&0xff; b[(*i)++]=(v>>8)&0xff; }
static void pF32(unsigned char* b, size_t* i, float f){ unsigned u; memcpy(&u,&f,4); pU32(b,i,u); }

/* section header: chromId, chromStart, chromEnd, step, span, type, res, itemCount */
static void hdr(unsigned char* b, size_t* i, unsigned cId, unsigned cStart, unsigned cEnd,
                unsigned step, unsigned span, unsigned type, unsigned items){
  pU32(b,i,cId); pU32(b,i,cStart); pU32(b,i,cEnd); pU32(b,i,step); pU32(b,i,span);
  b[(*i)++]=type; b[(*i)++]=0; pU16(b,i,items);
}

/* collect emitted intervals */
#define MAXC 64
static long cs[MAXC], ce[MAXC]; static float cv[MAXC]; static int cn;
static void collect(long s, long e, float v, void* ud){ (void)ud; cs[cn]=s; ce[cn]=e; cv[cn]=v; cn++; }

static int fails = 0;
static void expect(const char* name, int got, int want){
  if (got != want){ printf("FAIL %s: got %d intervals, want %d\n", name, got, want); fails++; }
}

int main(void){
  unsigned char buf[512];
  size_t n;

  /* --- fixedStep: cId=1 start=100 step=10 span=5, values 1,2,3 -> intervals
     [100,105)=1 [110,115)=2 [120,125)=3.  Query [108,122) overlaps last two. */
  n=0; hdr(buf,&n,1,100,130,10,5,SECT_FIXEDSTEP,3);
  pF32(buf,&n,1.0f); pF32(buf,&n,2.0f); pF32(buf,&n,3.0f);
  cn=0; emitSection(buf,n,1,108,122,collect,NULL);
  expect("fixedStep count",cn,2);
  if(cn==2){ assert(cs[0]==110&&ce[0]==115&&cv[0]==2.0f); assert(cs[1]==120&&ce[1]==125&&cv[1]==3.0f); }

  /* --- varStep: span=8, items (start=200,v=4) (start=250,v=5) -> [200,208) [250,258).
     Query [205,255) overlaps both. */
  n=0; hdr(buf,&n,1,200,258,0,8,SECT_VARSTEP,2);
  pU32(buf,&n,200); pF32(buf,&n,4.0f); pU32(buf,&n,250); pF32(buf,&n,5.0f);
  cn=0; emitSection(buf,n,1,205,255,collect,NULL);
  expect("varStep count",cn,2);
  if(cn==2){ assert(cs[0]==200&&ce[0]==208&&cv[0]==4.0f); assert(cs[1]==250&&ce[1]==258&&cv[1]==5.0f); }

  /* --- bedGraph: (300,320,6) (400,410,7); query [318,405) overlaps both, and the
     half-open edge: query end at exactly a start / query start at exactly an end. */
  n=0; hdr(buf,&n,1,300,410,0,0,SECT_BEDGRAPH,2);
  pU32(buf,&n,300); pU32(buf,&n,320); pF32(buf,&n,6.0f);
  pU32(buf,&n,400); pU32(buf,&n,410); pF32(buf,&n,7.0f);
  cn=0; emitSection(buf,n,1,318,405,collect,NULL);
  expect("bedGraph count",cn,2);

  /* half-open edges. A query sitting exactly in the gap [320,400) touches the end
     of [300,320) and the start of [400,410) but overlaps NEITHER. A one-base
     window at an item's last base overlaps; at an item's start-1 it does not. */
  cn=0; emitSection(buf,n,1,320,400,collect,NULL); expect("gap touches both edges -> none",cn,0);
  cn=0; emitSection(buf,n,1,319,320,collect,NULL); expect("last base of [300,320)",cn,1);
  if(cn==1) assert(cv[0]==6.0f);
  cn=0; emitSection(buf,n,1,400,401,collect,NULL); expect("first base of [400,410)",cn,1);
  if(cn==1) assert(cv[0]==7.0f);
  cn=0; emitSection(buf,n,1,290,300,collect,NULL); expect("just before [300,320) -> none",cn,0);

  /* wrong chrom id -> nothing */
  cn=0; emitSection(buf,n,2,300,410,collect,NULL); expect("wrong chrom",cn,0);

  /* unknown section type -> decode error (-1) */
  n=0; hdr(buf,&n,1,0,10,0,5,99,1); pF32(buf,&n,1.0f);
  if (emitSection(buf,n,1,0,10,collect,NULL) != -1){ printf("FAIL unknown-type not rejected\n"); fails++; }

  printf(fails ? "%d FAILED\n" : "section_unit: ALL PASS\n", fails);
  return fails ? 1 : 0;
}
