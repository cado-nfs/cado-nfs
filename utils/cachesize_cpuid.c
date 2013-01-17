#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "macros.h"

#ifdef HAVE_GCC_STYLE_AMD64_ASM

//#define DEBUG_CACHESIZE_CPUID

static const int UNKNOWN=1;
static const int INTEL=1;
static const int AMD=2;

static const uint32_t mask0 = 0xFF;
static const uint32_t mask1 = 0xFF00;
static const uint32_t mask2 = 0xFF0000;
static const uint32_t mask3 = 0xFF000000;
static const uint32_t mask_31_22 = 0xFFC00000;
static const uint32_t mask_21_12 = 0x3FF000;
static const uint32_t mask_11_0 = 0xFFF;

static inline uint32_t byte0(uint32_t x) {
  return (x & mask0);
}
static inline uint32_t byte1(uint32_t x) {
  return (x & mask1)>>8;
}
static inline uint32_t byte2(uint32_t x) {
  return (x & mask2)>>16;
}
static inline uint32_t byte3(uint32_t x) {
  return (x & mask3)>>24;
}
static inline uint32_t bits_31_22(uint32_t x) {
  return (x & mask_31_22)>>22;
}
static inline uint32_t bits_21_12(uint32_t x) {
  return (x & mask_21_12)>>12;
}
static inline uint32_t bits_11_0(uint32_t x) {
  return (x & mask_11_0);
}


void cpuid(uint32_t res[4], uint32_t op) {
    volatile uint64_t saved[4];
    volatile uint32_t pre_res[4];
    __asm__ __volatile__(
            "movq %%rax, 0(%[src])\n"
            "movq %%rbx, 8(%[src])\n"
            "movq %%rcx, 16(%[src])\n"
            "movq %%rdx, 24(%[src])\n"
            "movl %[op], %%eax\n"
            "cpuid\n"
            "movl %%eax, 0(%[res])\n"
            "movl %%ebx, 4(%[res])\n"
            "movl %%ecx, 8(%[res])\n"
            "movl %%edx, 12(%[res])\n"
            "movq 0(%[src]), %%rax\n"
            "movq 8(%[src]), %%rbx\n"
            "movq 16(%[src]), %%rcx\n"
            "movq 24(%[src]), %%rdx\n"
            :
            : [op]"m"(op), [res]"D"(pre_res), [src]"S"(saved)
            : "memory"
            );
    res[0] = pre_res[0];
    res[1] = pre_res[1];
    res[2] = pre_res[2];
    res[3] = pre_res[3];
}


// str should be 13 byte long, at least.
void vendor(char *str) {
  uint32_t abcd[4], x;
  cpuid(abcd, 0);
  x = abcd[1];
  str[0] = x & 255; x >>= 8;
  str[1] = x & 255; x >>= 8;
  str[2] = x & 255; x >>= 8;
  str[3] = x & 255; 
  x = abcd[3];
  str[4] = x & 255; x >>= 8;
  str[5] = x & 255; x >>= 8;
  str[6] = x & 255; x >>= 8;
  str[7] = x & 255; 
  x = abcd[2];
  str[8] = x & 255; x >>= 8;
  str[9] = x & 255; x >>= 8;
  str[10] = x & 255; x >>= 8;
  str[11] = x & 255; 
  str[12] = '\0';
}

int brand() {
  char *str = malloc(sizeof(char) * 13);
  vendor(str);
  if (strcmp(str, "AuthenticAMD")==0) {
    free(str);
    return AMD;
  }
  if (strcmp(str, "GenuineIntel")==0) {
    free(str);
    return INTEL;
  }
  free(str);
  return UNKNOWN;
}

typedef struct {
  int L1Data_size;  		// in KB
  int L1Data_assoc;
  int L1Data_line;  		// in B
  int DataTLB_pagesize; 	// in KB
  int DataTLB_entries;
  int DataTLB_assoc;
  int L1Instr_size;  		// in KB
  int L1Instr_assoc;
  int L1Instr_line;  		// in B
  int InstrTLB_pagesize; 	// in KB
  int InstrTLB_entries;
  int InstrTLB_assoc;
  int L2_size; 			// in KB
  int L3_size; 			// in KB
} cache_data_t;

void init_cache_data(cache_data_t *data) {
  data->L1Data_size = -1;  	
  data->L1Data_assoc = -1;
  data->L1Data_line = -1;  	
  data->DataTLB_pagesize = -1; 
  data->DataTLB_entries = -1;
  data->DataTLB_assoc = -1;
  data->L1Instr_size = -1;  	
  data->L1Instr_assoc = -1;
  data->L1Instr_line = -1;  	
  data->InstrTLB_pagesize = -1;
  data->InstrTLB_entries = -1;
  data->InstrTLB_assoc = -1;
  data->L2_size = -1; 		
  data->L3_size = -1; 		
}

void print_cache_data(cache_data_t *data) {
  if (data->L1Data_size != -1) printf("L1Data_size (KB) = %d\n", data->L1Data_size);  	
  if (data->L1Data_assoc != -1) printf("L1Data_assoc = %d\n", data->L1Data_assoc);
  if (data->L1Data_line != -1) printf("L1Data_line (B) = %d\n", data->L1Data_line);  	
  if (data->DataTLB_pagesize != -1) printf("DataTLB_pagesize (KB) = %d\n", data->DataTLB_pagesize); 
  if (data->DataTLB_entries != -1) printf("DataTLB_entries = %d\n", data->DataTLB_entries);
  if (data->DataTLB_assoc != -1) printf("DataTLB_assoc = %d\n", data->DataTLB_assoc);
  if (data->L1Instr_size != -1) printf("L1Instr_size (KB) = %d\n", data->L1Instr_size);  	
  if (data->L1Instr_assoc != -1) printf("L1Instr_assoc = %d\n", data->L1Instr_assoc);
  if (data->L1Instr_line != -1) printf("L1Instr_line (B) = %d\n", data->L1Instr_line);  	
  if (data->InstrTLB_pagesize != -1) printf("InstrTLB_pagesize (KB) = %d\n", data->InstrTLB_pagesize);
  if (data->InstrTLB_entries != -1) printf("InstrTLB_entries = %d\n", data->InstrTLB_entries);
  if (data->InstrTLB_assoc != -1) printf("InstrTLB_assoc = %d\n", data->InstrTLB_assoc);
  if (data->L2_size != -1) printf("L2_size (KB) = %d\n", data->L2_size); 		
  if (data->L3_size != -1) printf("L3_size (KB) = %d\n", data->L3_size); 		
}


// Taken from Table 3-25 in Intel ref manual number 253666
// pages: 3-198 and following, in particular Table 3-22.
void
update_intel_byte (uint32_t c, cache_data_t *data) {
  if (c == 0)
    return;
  switch (c) {
   case 1:
    data->InstrTLB_pagesize = 4;
    data->InstrTLB_entries = 32;
    data->InstrTLB_assoc = 4;
    break;
   case 2:
    // for 4M TLB
    break;
   case 3:
    data->DataTLB_pagesize = 4;
    data->DataTLB_entries = 64;
    data->DataTLB_assoc = 4;
    break;
   case 4:
   case 5:
    // for 4M TLB
    break;
   case 6:
    data->L1Instr_size = 8;
    data->L1Instr_assoc = 4;
    data->L1Instr_line = 32;
    break;
   case 8:
    data->L1Instr_size = 16;
    data->L1Instr_assoc = 4;
    data->L1Instr_line = 32;
    break;
   case 0xA:
    data->L1Data_size = 8;
    data->L1Data_assoc = 2;
    data->L1Data_line = 32;
    break;
   case 0xB:
    // for 4M TLB
    break;
   case 0xC:
    data->L1Data_size = 16;
    data->L1Data_assoc = 4;
    data->L1Data_line = 32;
    break;
   case 0x22:
    data->L3_size = 512;
    break;
   case 0x23:
    data->L3_size = 1024;
    break;
   case 0x25:
    data->L3_size = 2048;
    break;
   case 0x29:
    data->L3_size = 4096;
    break;
   case 0x2C:
    data->L1Data_size = 32;
    data->L1Data_assoc = 8;
    data->L1Data_line = 64;
    break;
   case 0x30:
    data->L1Instr_size = 32;
    data->L1Instr_assoc = 8;
    data->L1Instr_line = 64;
    break;
   case 0x41:
    data->L2_size = 128;
    break;
   case 0x42:
    data->L2_size = 256;
    break;
   case 0x43:
    data->L2_size = 512;
    break;
   case 0x44:
    data->L2_size = 1024;
    break;
   case 0x45:
    data->L2_size = 2048;
    break;
   case 0x46:
    data->L3_size = 4096;
    break;
   case 0x47:
    data->L3_size = 8192;
    break;
   case 0x48:
    data->L2_size = 3072;
    break;
   case 0x49:
    data->L2_size = 4096;
    break;
   case 0x4A:
    data->L3_size = 6144;
    break;
   case 0x4B:
    data->L3_size = 8192;
    break;
   case 0x4C:
    data->L3_size = 12288;
    break;
   case 0x4D:
    data->L3_size = 16384;
    break;
   case 0x4E:
    data->L2_size = 6144;
    break;
   case 0x50:
    data->InstrTLB_pagesize = 4;
    data->InstrTLB_entries = 64;
    break;
   case 0x51:
    data->InstrTLB_pagesize = 4;
    data->InstrTLB_entries = 128;
    break;
   case 0x52:
    data->InstrTLB_pagesize = 4;
    data->InstrTLB_entries = 256;
    break;
   case 0x56:
   case 0x57:
    // for 4MB TLB
    break;
   case 0x5B:
    data->DataTLB_pagesize = 4;
    data->DataTLB_entries = 64;
    break;
   case 0x5C:
    data->DataTLB_pagesize = 4;
    data->DataTLB_entries = 128;
    break;
   case 0x5D:
    data->DataTLB_pagesize = 4;
    data->DataTLB_entries = 256;
    break;
   case 0x60:
    data->L1Data_size = 16;
    data->L1Data_assoc = 8;
    data->L1Data_line = 64;
    break;
   case 0x66:
    data->L1Data_size = 8;
    data->L1Data_assoc = 4;
    data->L1Data_line = 64;
    break;
   case 0x67:
    data->L1Data_size = 16;
    data->L1Data_assoc = 4;
    data->L1Data_line = 64;
    break;
   case 0x68:
    data->L1Data_size = 32;
    data->L1Data_assoc = 4;
    data->L1Data_line = 64;
    break;
   case 0x70:
   case 0x71:
   case 0x72:
    // trace cache
    break;
   case 0x78:
    data->L2_size = 1024;
    break;
   case 0x79:
    data->L2_size = 128;
    break;
   case 0x7A:
    data->L2_size = 256;
    break;
   case 0x7B:
    data->L2_size = 512;
    break;
   case 0x7C:
    data->L2_size = 1024;
    break;
   case 0x7D:
    data->L2_size = 2048;
    break;
   case 0x7F:
    data->L2_size = 512;
    break;
   case 0x82:
    data->L2_size = 256;
    break;
   case 0x83:
    data->L2_size = 512;
    break;
   case 0x84:
    data->L2_size = 1024;
    break;
   case 0x85:
    data->L2_size = 2048;
    break;
   case 0x86:
    data->L2_size = 512;
    break;
   case 0x87:
    data->L2_size = 1024;
    break;
   case 0xB0:
    data->InstrTLB_pagesize = 4;
    data->InstrTLB_entries = 128;
    data->InstrTLB_assoc = 4;
    break;
   case 0xB1:
    // for TLB of 4MB
    break;
   case 0xB3:
    data->DataTLB_pagesize = 4;
    data->DataTLB_entries = 128;
    data->DataTLB_assoc = 4;
    break;
   case 0xB4:
    data->DataTLB_pagesize = 4;
    data->DataTLB_entries = 256;
    data->DataTLB_assoc = 4;
    break;
   case 0xF0:
   case 0xF1:
    // for prefetching
    break;
   case 0xFF:
   {
     // cpuid 4
     uint32_t res[4];
     cpuid(res, 4);
#ifdef DEBUG_CACHESIZE_CPUID
     int j = 0;
     for (j = 0; j < 4; j ++)
       printf ("res[%d]: %08X\n", j, res[j]);
     exit(1);
#endif
     // (EBX[31:22] + 1) * (EBX[21:12] + 1) * (EBX[11:0] + 1) * (ECX + 1)
     data->L1Data_size = (bits_31_22(res[1]) + 1) *
       (bits_21_12(res[1]) + 1) * (bits_11_0(res[1]) + 1) * (res[2] + 1) / 1024;
     break;
   }
   default:
    break;
  }
}


int
print_intel_cache (int verbose)
{
  uint32_t res[4];

  cpuid(res, 0);
  if (res[0] < 2) {
    if (verbose)
      printf("No cache data available\n");
    return -1;
  }
  cpuid(res, 2);
  if (byte0(res[0]) != 1) {
    if (verbose)
      printf("Need multiple calls to cpuid(2). Sorry, not implemented!\n");
    return -1;
  }
  cache_data_t data;
  init_cache_data(&data);
  int i;
  for (i = 0; i < 4; ++i) {
    if (res[i] & 0x80000000)
      continue;
    update_intel_byte(byte0(res[i]), &data);
    update_intel_byte(byte1(res[i]), &data);
    update_intel_byte(byte2(res[i]), &data);
    update_intel_byte(byte3(res[i]), &data);
  }

  if (verbose)
    print_cache_data (&data);

  return ((&data)->L1Data_size);  	
}


int print_amd_cache(int verbose) {
  uint32_t res[4];
  cpuid(res, 0x80000000);
  if (res[0] < 0x80000005) {
    if (verbose)
      printf("No cache data available\n");
    return -1;
  }
  cache_data_t data;
  init_cache_data(&data);

  cpuid(res, 0x80000005);
  data.L1Data_size = byte3(res[2]);
  data.L1Data_line = byte0(res[2]);
  data.L1Data_assoc = byte2(res[2]);
  data.DataTLB_pagesize = 4; // always the same for AMD, it seems.
  data.DataTLB_entries = byte2(res[1]);
  data.DataTLB_assoc = byte3(res[1]);
  data.L1Instr_size = byte3(res[3]);
  data.L1Instr_line = byte0(res[3]);
  data.L1Instr_assoc = byte2(res[3]);
  data.InstrTLB_pagesize = 4; // always the same for AMD, it seems.
  data.InstrTLB_entries = byte0(res[1]);
  data.InstrTLB_assoc = byte1(res[1]);

  cpuid(res, 0x80000000);
  if (res[0] >= 0x80000006) {
    cpuid(res, 0x80000006);
    data.L2_size = res[2]>>16;
    data.L3_size = 512*(res[3]>>18);
  }
  
  if (verbose)
    print_cache_data(&data);
  
  return ((&data)->L1Data_size);
}


/* return -1 if cpuid() failed, you need to seek for alternative
   ways to probe it */
int
cachesize_cpuid (int verbose)
{

  int brd, ret = -1;

  brd = brand ();
  if (brd == AMD)
    ret = print_amd_cache (verbose);
  else if (brd == INTEL)
    ret = print_intel_cache (verbose);
  else
    {
      if (verbose) {
        fprintf (stderr, "Warn, unknown architecture in cachesize_cpuid()\n");
      }
      return -1; // continue here.
    }
  return ret * 1024;
}

#else

int
cachesize_cpuid (int verbose MAYBE_UNUSED)
{
  return -1;
}

#endif
