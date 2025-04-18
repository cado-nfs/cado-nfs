#ifndef CADO_VERSION_INFO_H
#define CADO_VERSION_INFO_H

#ifdef  __cplusplus
extern "C" {
#endif

/* These symbols are provided by the cmake-generated file version_info.c,
 * from version_info.c.in
 */
extern const char * cado_revision_string;

extern unsigned int cado_version_major;
extern unsigned int cado_version_minor;
extern const char * cado_version_string;

/* This symbol is provided by the script-generated file modified_files.c
 */
extern const char * cado_modified_files;

#ifdef  __cplusplus
}
#endif

#endif	/* CADO_VERSION_INFO_H */
