#ifndef PORTOPT_PORTOPT_H
#define PORTOPT_PORTOPT_H

#include <stddef.h>
#include <stdint.h>

#if defined(_WIN32)
#if defined(PO_BUILDING)
#define PO_API __declspec(dllexport)
#else
#define PO_API __declspec(dllimport)
#endif
#else
#define PO_API __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C"
{
#endif

#define PO_ABI_VERSION 0

PO_API int32_t po_abi_version(void);

/* Library version, "major.minor.patch". The pointer is valid for the lifetime of the process. */
PO_API const char* po_version(void);

typedef int32_t po_status;
#define PO_OK 0
#define PO_INVALID_ARGUMENT 1
#define PO_INVALID_HANDLE 2
#define PO_UNSUPPORTED 3
#define PO_OUT_OF_MEMORY 4
#define PO_INTERNAL_ERROR 5

/* Message of the last failed call on the calling thread ("" after a success). The pointer stays valid
   until the next PortOpt call on that thread. */
PO_API const char* po_last_error(void);

#ifdef __cplusplus
}
#endif

#endif
