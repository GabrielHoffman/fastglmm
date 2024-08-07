
// Adapted from RhpcBLASctl R package on CRAN
// Aug 7, 2024
// Gabriel Hoffman
// Icahn School of Medicine at Mount Sinai

#define OPTION_WIN_RBLAS_NAME "RhpcBLASctl.win.Rblas.name"
#define OPTION_WIN_RBLAS_TYPE "RhpcBLASctl.win.Rblas.type"

// static char *getRblasName(void)
// {
//   static char Rblas[1024];
//   SEXP ret = GetOption1(install(OPTION_WIN_RBLAS_NAME));
//   if(ret==R_NilValue||XLENGTH(ret)<1||TYPEOF(ret)!=STRSXP) return("Rblas.dll");
//   strncpy(Rblas, CHAR(STRING_ELT(ret,0)), sizeof(Rblas)-1);
//   return(Rblas);
// }


#ifdef WIN32
#include <windows.h>
static  void* DLOPEN(void)
{
  return((void*)LoadLibrary(getRblasName()));
}
static  FARPROC DLSYM(void *handle, const char *symbol)
{
  FARPROC proc = GetProcAddress((HMODULE)handle,(LPCSTR)symbol);
  if(proc != NULL)
    Rprintf("detected function %s\n", symbol);
  return(proc);
}
#define EXPDLSYM(_res,_cast,_handle,_symbol) (_res = (_cast)DLSYM(_handle,_symbol)) 
static  int DLCLOSE(void *handle)
{
  return(FreeLibrary((HMODULE)handle));
}
#else
#include <dlfcn.h>
static  void* DLOPEN(void)
{
  return(dlopen(NULL, RTLD_NOW|RTLD_GLOBAL));
}
static  void* DLSYM(void *handle, const char *symbol)
{
  void* proc;
  *(void**)(&proc)= dlsym(handle,symbol);
  return(proc);
}
#define EXPDLSYM(_res,_cast,_handle,_symbol) (*(void**)(&_res) =DLSYM(_handle,_symbol)) 
static  int DLCLOSE(void *handle)
{
  return(dlclose(handle));
}
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef WIN32
#include "config.h"
#endif

#include <R.h>
#include <Rinternals.h>



#ifdef HAVE_SYS_SYSINFO_H
#include <sys/sysinfo.h>
#endif

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif

#ifdef HAVE_SYS_SYSCTL_H
#if !defined(__linux__)
#include <sys/sysctl.h>
#endif /* !defined(__linux__) */
#endif

#if defined(WIN32)
#if !defined(LTP_PC_SMT)
typedef enum _PROCESSOR_CACHE_TYPE {
  CacheUnified,CacheInstruction,CacheData,CacheTrace
} PROCESSOR_CACHE_TYPE;
typedef struct _CACHE_DESCRIPTOR {
  BYTE Level;
  BYTE Associativity;
  WORD LineSize;
  DWORD Size;
  PROCESSOR_CACHE_TYPE Type;
} CACHE_DESCRIPTOR,*PCACHE_DESCRIPTOR;
#endif
#if !defined(CACHE_FULLY_ASSOCIATIVE)
typedef enum _LOGICAL_PROCESSOR_RELATIONSHIP {
  RelationProcessorCore,RelationNumaNode,RelationCache
} LOGICAL_PROCESSOR_RELATIONSHIP;
typedef struct _SYSTEM_LOGICAL_PROCESSOR_INFORMATION {
  ULONG_PTR ProcessorMask;
  LOGICAL_PROCESSOR_RELATIONSHIP Relationship;
  union {
    struct {
      BYTE Flags;
    } ProcessorCore;
    struct {
      DWORD NodeNumber;
    } NumaNode;
    CACHE_DESCRIPTOR Cache;
    ULONGLONG Reserved[2];
  };
} SYSTEM_LOGICAL_PROCESSOR_INFORMATION,*PSYSTEM_LOGICAL_PROCESSOR_INFORMATION;
#endif
#endif

#define MKL_ALL  0
#define MKL_BLAS 1
#define MKL_FFT  2
#define MKL_VML  3

#define    STR_GOTO_GET_NUM_PROCS         "goto_get_num_procs"
#define    STR_GOTO_SET_NUM_THREADS       "goto_set_num_threads"
#define    STR_MKL_DOMAIN_GET_MAX_THREADS "mkl_domain_get_max_threads"
#define    STR_MKL_DOMAIN_SET_NUM_THREADS "mkl_domain_set_num_threads"
#define    STR_MKL_GET_MAX_THREADS        "mkl_get_max_threads"
#define    STR_MKL_SET_NUM_THREADS        "mkl_set_num_threads"
#define    STR_ACMLGETMAXTHREADS          "acmlgetmaxthreads"
#define    STR_ACMLSETNUMTHREADS          "acmlsetnumthreads"
#define    STR_BLI_THREAD_GET_NUM_THREADS "bli_thread_get_num_threads"
#define    STR_BLI_THREAD_SET_NUM_THREADS "bli_thread_set_num_threads"

typedef int  (*GOTO_GET_NUM_PROCS)        (void);
typedef void (*GOTO_SET_NUM_THREADS)      (int);
typedef int  (*MKL_DOMAIN_GET_MAX_THREADS)(int*);
typedef int  (*MKL_DOMAIN_SET_NUM_THREADS)(int*,int*);
typedef int  (*MKL_GET_MAX_THREADS)       (void);
typedef void (*MKL_SET_NUM_THREADS)       (int*);
typedef int  (*ACMLGETMAXTHREADS)         (void);
typedef void (*ACMLSETNUMTHREADS)         (int);
typedef int  (*BLI_THREAD_GET_NUM_THREADS)(void);
typedef void (*BLI_THREAD_SET_NUM_THREADS)(int);
 
static         GOTO_GET_NUM_PROCS          goto_get_num_procs          = NULL;
static         GOTO_SET_NUM_THREADS        goto_set_num_threads        = NULL;
static         MKL_DOMAIN_GET_MAX_THREADS  mkl_domain_get_max_threads  = NULL;
static         MKL_DOMAIN_SET_NUM_THREADS  mkl_domain_set_num_threads  = NULL;
static         MKL_GET_MAX_THREADS         mkl_get_max_threads         = NULL;
static         MKL_SET_NUM_THREADS         mkl_set_num_threads         = NULL;
static         ACMLGETMAXTHREADS           acmlgetmaxthreads           = NULL;
static         ACMLSETNUMTHREADS           acmlsetnumthreads           = NULL;
static         BLI_THREAD_GET_NUM_THREADS  bli_thread_get_num_threads  = NULL;
static         BLI_THREAD_SET_NUM_THREADS  bli_thread_set_num_threads  = NULL;


#ifdef WIN32
typedef BOOL (WINAPI *GLPI)(PSYSTEM_LOGICAL_PROCESSOR_INFORMATION, PDWORD);
#endif



int blas_set_num_threads(int threads){

  void *dlh = DLOPEN();

  if(dlh==NULL){
    // error("Failed to acquire BLAS handle.");
    return -1;
  }

  if (threads < 1) threads=1; /* minimum 1 */
  
  if(      NULL != ( EXPDLSYM(goto_set_num_threads,       GOTO_SET_NUM_THREADS,       dlh, STR_GOTO_SET_NUM_THREADS))){
    goto_set_num_threads(threads);
  }
  else if( NULL != ( EXPDLSYM(mkl_domain_set_num_threads, MKL_DOMAIN_SET_NUM_THREADS, dlh, STR_MKL_DOMAIN_SET_NUM_THREADS))){
    int   type = MKL_BLAS;
    mkl_domain_set_num_threads(&threads,&type);
  }
  else if(NULL != ( EXPDLSYM(mkl_set_num_threads,         MKL_SET_NUM_THREADS,        dlh, STR_MKL_SET_NUM_THREADS))){
    mkl_set_num_threads(&threads);
  }
  else if( NULL != ( EXPDLSYM(acmlsetnumthreads,          ACMLSETNUMTHREADS,          dlh, STR_ACMLSETNUMTHREADS))){
    acmlsetnumthreads(threads);
  }
  else if( NULL != ( EXPDLSYM(bli_thread_set_num_threads, BLI_THREAD_SET_NUM_THREADS, dlh, STR_BLI_THREAD_SET_NUM_THREADS))){
    bli_thread_set_num_threads(threads);
  }

  DLCLOSE(dlh);
  return -1;
}


