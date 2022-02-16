#pragma unknown_control_flow(_setjmp)
#pragma unknown_control_flow(__setjmp)
#pragma unknown_control_flow(__sigsetjmp)
#pragma unknown_control_flow(vfork)

/***************************************************************
 * BIT UTILITY for C and FORTRAN    
 * 2003.10.24 recieve from SATO-san, SATELLITE-GROUP, NWP, JMA.
 ***************************************************************/

/* Define ISO C stdio on top of C++ iostreams.
   Copyright (C) 1991,1994-2004,2005,2006,2007 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */

/*
 *	ISO C99 Standard: 7.19 Input/output	<stdio.h>
 */


/* features.h Copyright FUJITSU LIMITED 2008-2010 */

/* common definition to all modes */




/* Copyright (C) 1991,1992,1993,1995-2006,2007 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */


/* These are defined by the user (or the compiler)
   to specify the desired environment:

   __STRICT_ANSI__	ISO Standard C.
   _ISOC99_SOURCE	Extensions to ISO C89 from ISO C99.
   _POSIX_SOURCE	IEEE Std 1003.1.
   _POSIX_C_SOURCE	If ==1, like _POSIX_SOURCE; if >=2 add IEEE Std 1003.2;
			if >=199309L, add IEEE Std 1003.1b-1993;
			if >=199506L, add IEEE Std 1003.1c-1995;
			if >=200112L, all of IEEE 1003.1-2004
   _XOPEN_SOURCE	Includes POSIX and XPG things.  Set to 500 if
			Single Unix conformance is wanted, to 600 for the
			upcoming sixth revision.
   _XOPEN_SOURCE_EXTENDED XPG things and X/Open Unix extensions.
   _LARGEFILE_SOURCE	Some more functions for correct standard I/O.
   _LARGEFILE64_SOURCE	Additional functionality from LFS for large files.
   _FILE_OFFSET_BITS=N	Select default filesystem interface.
   _BSD_SOURCE		ISO C, POSIX, and 4.3BSD things.
   _SVID_SOURCE		ISO C, POSIX, and SVID things.
   _ATFILE_SOURCE	Additional *at interfaces.
   _GNU_SOURCE		All of the above, plus GNU extensions.
   _REENTRANT		Select additionally reentrant object.
   _THREAD_SAFE		Same as _REENTRANT, often used by other systems.
   _FORTIFY_SOURCE	If set to numeric value > 0 additional security
			measures are defined, according to level.

   The `-ansi' switch to the GNU C compiler defines __STRICT_ANSI__.
   If none of these are defined, the default is to have _SVID_SOURCE,
   _BSD_SOURCE, and _POSIX_SOURCE set to one and _POSIX_C_SOURCE set to
   200112L.  If more than one of these are defined, they accumulate.
   For example __STRICT_ANSI__, _POSIX_SOURCE and _POSIX_C_SOURCE
   together give you ISO C, 1003.1, and 1003.2, but nothing else.

   These are defined by this file and are used by the
   header files to decide what to declare or define:

   __USE_ISOC99		Define ISO C99 things.
   __USE_ISOC95		Define ISO C90 AMD1 (C95) things.
   __USE_POSIX		Define IEEE Std 1003.1 things.
   __USE_POSIX2		Define IEEE Std 1003.2 things.
   __USE_POSIX199309	Define IEEE Std 1003.1, and .1b things.
   __USE_POSIX199506	Define IEEE Std 1003.1, .1b, .1c and .1i things.
   __USE_XOPEN		Define XPG things.
   __USE_XOPEN_EXTENDED	Define X/Open Unix things.
   __USE_UNIX98		Define Single Unix V2 things.
   __USE_XOPEN2K        Define XPG6 things.
   __USE_LARGEFILE	Define correct standard I/O things.
   __USE_LARGEFILE64	Define LFS things with separate names.
   __USE_FILE_OFFSET64	Define 64bit interface as default.
   __USE_BSD		Define 4.3BSD things.
   __USE_SVID		Define SVID things.
   __USE_MISC		Define things common to BSD and System V Unix.
   __USE_ATFILE		Define *at interfaces and AT_* constants for them.
   __USE_GNU		Define GNU extensions.
   __USE_REENTRANT	Define reentrant/thread-safe *_r functions.
   __USE_FORTIFY_LEVEL	Additional security measures used, according to level.
   __FAVOR_BSD		Favor 4.3BSD things in cases of conflict.

   The macros `__GNU_LIBRARY__', `__GLIBC__', and `__GLIBC_MINOR__' are
   defined by this file unconditionally.  `__GNU_LIBRARY__' is provided
   only for compatibility.  All new code should use the other symbols
   to test for features.

   All macros listed above as possibly being defined by this file are
   explicitly undefined if they are not explicitly defined.
   Feature-test macros that are not defined by the user or compiler
   but are implied by the other feature-test macros defined (or by the
   lack of any definitions) are defined by the file.  */


/* Undefine everything, so we get a clean slate.  */

/* Suppress kernel-name space pollution unless user expressedly asks
   for it.  */

/* Always use ISO C things.  */

/* Convenience macros to test the versions of glibc and gcc.
   Use them like this:
   #if __GNUC_PREREQ (2,8)
   ... code requiring gcc 2.8 or later ...
   #endif
   Note - they won't work for gcc1 or glibc1, since the _MINOR macros
   were not defined then.  */


/* If _BSD_SOURCE was defined by the user, favor BSD over POSIX.  */

/* If _GNU_SOURCE was defined by the user, turn on all the other features.  */

/* If nothing (other than _GNU_SOURCE) is defined,
   define _BSD_SOURCE and _SVID_SOURCE.  */

/* This is to enable the ISO C99 extension.  Also recognize the old macro
   which was used prior to the standard acceptance.  This macro will
   eventually go away and the features enabled by default once the ISO C99
   standard is widely adopted.  */

/* This is to enable the ISO C90 Amendment 1:1995 extension.  */

/* If none of the ANSI/POSIX macros are defined, use POSIX.1 and POSIX.2
   (and IEEE Std 1003.1b-1993 unless _XOPEN_SOURCE is defined).  */

















/* We do support the IEC 559 math functionality, real and complex.  */

/* wchar_t uses ISO 10646-1 (2nd ed., published 2000-09-15) / Unicode 3.1.  */

/* This macro indicates that the installed library is the GNU C Library.
   For historic reasons the value now is 6 and this will stay from now
   on.  The use of this variable is deprecated.  Use __GLIBC__ and
   __GLIBC_MINOR__ now (see below) when you want to test for a specific
   GNU C library version and use the values in <gnu/lib-names.h> to get
   the sonames of the shared libraries.  */

/* Major and minor version number of the GNU C library package.  Use
   these macros to test for features in specific releases.  */


/* Decide whether a compiler supports the long long datatypes.  */

/* This is here only because every header file already includes this one.  */
/* Copyright (C) 1992-2001, 2002, 2004, 2005, 2006, 2007
   Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */


/* We are almost always included from features.h. */

/* The GNU libc does not support any K&R compilers or the traditional mode
   of ISO C compilers anymore.  Check for some of the combinations not
   anymore supported.  */

/* Some user header file might have defined this before.  */






/* These two macros are not used in glibc anymore.  They are kept here
   only because some other projects expect the macros to be defined.  */

/* For these things, GCC behaves the ANSI way normally,
   and the non-ANSI way under -traditional.  */


/* This is not a typedef so `const __ptr_t' does the right thing.  */


/* C++ needs to know that types and declarations are C, not C++.  */


/* The standard library needs the functions from the ISO C90 standard
   in the std namespace.  At the same time we want to be safe for
   future changes and we include the ISO C99 code in the non-standard
   namespace __c99.  The C++ wrapper header take case of adding the
   definitions to the global namespace.  */
/* For compatibility we do not add the declarations into any
   namespace.  They will end up in the global namespace which is what
   old code expects.  */


/* Support for bounded pointers.  */


/* Fortify support.  */


/* Support for flexible arrays.  */


/* __asm__ ("xyz") is used throughout the headers to rename functions
   at the assembly language level.  This is wrapped by the __REDIRECT
   macro, in order to support compilers that can do this some other
   way.  When compilers don't support asm-names at all, we have to do
   preprocessor tricks instead (which don't have exactly the right
   semantics, but it's the best we can do).

   Example:
   int __REDIRECT(setpgrp, (__pid_t pid, __pid_t pgrp), setpgid); */


/* GCC has various useful declarations that can be made with the
   `__attribute__' syntax.  All of the ways we use this do fine if
   they are omitted for compilers that don't understand it. */

/* At some point during the gcc 2.96 development the `malloc' attribute
   for functions was introduced.  We don't want to use it unconditionally
   (although this would be possible) since it generates warnings.  */

/* At some point during the gcc 2.96 development the `pure' attribute
   for functions was introduced.  We don't want to use it unconditionally
   (although this would be possible) since it generates warnings.  */

/* At some point during the gcc 3.1 development the `used' attribute
   for functions was introduced.  We don't want to use it unconditionally
   (although this would be possible) since it generates warnings.  */

/* gcc allows marking deprecated functions.  */

/* At some point during the gcc 2.8 development the `format_arg' attribute
   for functions was introduced.  We don't want to use it unconditionally
   (although this would be possible) since it generates warnings.
   If several `format_arg' attributes are given for the same function, in
   gcc-3.0 and older, all but the last one are ignored.  In newer gccs,
   all designated arguments are considered.  */

/* At some point during the gcc 2.97 development the `strfmon' format
   attribute for functions was introduced.  We don't want to use it
   unconditionally (although this would be possible) since it
   generates warnings.  */

/* The nonull function attribute allows to mark pointer parameters which
   must not be NULL.  */

/* If fortification mode, we warn about unused results of certain
   function calls which can lead to problems.  */

/* Forces a function to be always inlined.  */

/* GCC 4.3 and above with -std=c99 or -std=gnu99 implements ISO C99
   inline semantics, unless -fgnu89-inline is used.  */

/* GCC 4.3 and above allow passing all anonymous arguments of an
   __extern_always_inline function to some other vararg function.  */

/* It is possible to compile containing GCC extensions even if GCC is
   run in pedantic mode if the uses are carefully marked using the
   `__extension__' keyword.  But this is not generally available before
   version 2.8.  */

/* __restrict is known in EGCS 1.2 and above. */

/* ISO C99 also allows to declare arrays as non-overlapping.  The syntax is
     array_name[restrict]
   GCC 3.1 supports this.  */

/* Determine the wordsize from the preprocessor defines.  */






/* If we don't have __REDIRECT, prototypes will be missing if
   __USE_FILE_OFFSET64 but not __USE_LARGEFILE[64]. */


/* Decide whether we can define 'extern inline' functions in headers.  */


/* This is here only because every header file already includes this one.
   Get the definitions of all the appropriate `__stub_FUNCTION' symbols.
   <gnu/stubs.h> contains `#define __stub_FUNCTION' when FUNCTION is a stub
   that will always return failure (and set errno to ENOSYS).  */
/* This file selects the right generated file of `__stub_FUNCTION' macros
   based on the architecture being compiled for.  */

/* Determine the wordsize from the preprocessor defines.  */




/* This file is automatically generated.
   It defines a symbol `__stub_FUNCTION' for each function
   in the C library which is a stub, meaning it will fail
   every time called, usually setting errno to ENOSYS.  */







/* stddef.h COPYRIGHT FUJITSU LIMITED 1998-2009                          */


 

 
 
typedef unsigned long size_t;
 





/* bits/types.h -- definitions of __*_t types underlying *_t types.
   Copyright (C) 2002, 2003, 2004, 2005, 2007 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */

/*
 * Never include this file directly; use <sys/types.h> instead.
 */


/* features.h Copyright FUJITSU LIMITED 2008-2010 */

/* Determine the wordsize from the preprocessor defines.  */




/* Convenience types.  */
typedef unsigned char __u_char;
typedef unsigned short int __u_short;
typedef unsigned int __u_int;
typedef unsigned long int __u_long;

/* Fixed-size types, underlying types depend on word size and compiler.  */
typedef signed char __int8_t;
typedef unsigned char __uint8_t;
typedef signed short int __int16_t;
typedef unsigned short int __uint16_t;
typedef signed int __int32_t;
typedef unsigned int __uint32_t;
typedef signed long int __int64_t;
typedef unsigned long int __uint64_t;

/* quad_t is also 64 bits.  */
typedef long int __quad_t;
typedef unsigned long int __u_quad_t;


/* The machine-dependent file <bits/typesizes.h> defines __*_T_TYPE
   macros for each of the OS types we define below.  The definitions
   of those macros must use the following macros for underlying types.
   We define __S<SIZE>_TYPE and __U<SIZE>_TYPE for the signed and unsigned
   variants of each of the following integer types on this machine.

	16		-- "natural" 16-bit type (always short)
	32		-- "natural" 32-bit type (always int)
	64		-- "natural" 64-bit type (long or long long)
	LONG32		-- 32-bit type, traditionally long
	QUAD		-- 64-bit type, always long long
	WORD		-- natural type of __WORDSIZE bits (int or long)
	LONGWORD	-- type of __WORDSIZE bits, traditionally long

   We distinguish WORD/LONGWORD, 32/LONG32, and 64/QUAD so that the
   conventional uses of `long' or `long long' type modifiers match the
   types we define, even when a less-adorned type would be the same size.
   This matters for (somewhat) portably writing printf/scanf formats for
   these types, where using the appropriate l or ll format modifiers can
   make the typedefs and the formats match up across all GNU platforms.  If
   we used `long' when it's 64 bits where `long long' is expected, then the
   compiler would warn about the formats not matching the argument types,
   and the programmer changing them to shut up the compiler would break the
   program's portability.

   Here we assume what is presently the case in all the GCC configurations
   we support: long long is always 64 bits, long is always word/address size,
   and int is always 32 bits.  */

/* No need to mark the typedef with __extension__.   */
/* bits/typesizes.h -- underlying types for *_t.  Linux/SPARC version.
   Copyright (C) 2002, 2003 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */



/* See <bits/types.h> for the meaning of these macros.  This file exists so
   that <bits/types.h> need not vary across different GNU platforms.  */


/* Number of descriptors that can fit in an `fd_set'.  */




typedef unsigned long int __dev_t;	/* Type of device numbers.  */
typedef unsigned int __uid_t;	/* Type of user identifications.  */
typedef unsigned int __gid_t;	/* Type of group identifications.  */
typedef unsigned long int __ino_t;	/* Type of file serial numbers.  */
typedef unsigned long int __ino64_t;	/* Type of file serial numbers (LFS).*/
typedef unsigned int __mode_t;	/* Type of file attribute bitmasks.  */
typedef unsigned int __nlink_t;	/* Type of file link counts.  */
typedef long int __off_t;	/* Type of file sizes and offsets.  */
typedef long int __off64_t;	/* Type of file sizes and offsets (LFS).  */
typedef int __pid_t;	/* Type of process identifications.  */
typedef struct { int __val[2]; } __fsid_t;	/* Type of file system IDs.  */
typedef long int __clock_t;	/* Type of CPU usage counts.  */
typedef unsigned long int __rlim_t;	/* Type for resource measurement.  */
typedef unsigned long int __rlim64_t;	/* Type for resource measurement (LFS).  */
typedef unsigned int __id_t;		/* General type for IDs.  */
typedef long int __time_t;	/* Seconds since the Epoch.  */
typedef unsigned int __useconds_t; /* Count of microseconds.  */
typedef int __suseconds_t; /* Signed count of microseconds.  */

typedef int __daddr_t;	/* The type of a disk address.  */
typedef long int __swblk_t;	/* Type of a swap block maybe?  */
typedef int __key_t;	/* Type of an IPC key.  */

/* Clock ID used in clock and timer functions.  */
typedef int __clockid_t;

/* Timer ID returned by `timer_create'.  */
typedef void * __timer_t;

/* Type to represent block size.  */
typedef long int __blksize_t;

/* Types from the Large File Support interface.  */

/* Type to count number of disk blocks.  */
typedef long int __blkcnt_t;
typedef long int __blkcnt64_t;

/* Type to count file system blocks.  */
typedef unsigned long int __fsblkcnt_t;
typedef unsigned long int __fsblkcnt64_t;

/* Type to count file system nodes.  */
typedef unsigned long int __fsfilcnt_t;
typedef unsigned long int __fsfilcnt64_t;

typedef long int __ssize_t; /* Type of a byte count, or error.  */

/* These few don't really vary by system, they always correspond
   to one of the other defined types.  */
typedef __off64_t __loff_t;	/* Type of file sizes and offsets (LFS).  */
typedef __quad_t *__qaddr_t;
typedef char *__caddr_t;

/* Duplicates info from stdint.h but this is used in unistd.h.  */
typedef long int __intptr_t;

/* Duplicate info from sys/socket.h.  */
typedef unsigned int __socklen_t;






/* Define outside of namespace so the C++ is happy.  */
struct _IO_FILE;


/* The opaque type of streams.  This is the definition used elsewhere.  */
typedef struct _IO_FILE FILE;






/* The opaque type of streams.  This is the definition used elsewhere.  */
typedef struct _IO_FILE __FILE;




/* Copyright (C) 1991-1995,1997-2006,2007 Free Software Foundation, Inc.
   This file is part of the GNU C Library.
   Written by Per Bothner <bothner@cygnus.com>.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.

   As a special exception, if you link the code in this file with
   files compiled with a GNU compiler to produce an executable,
   that does not cause the resulting executable to be covered by
   the GNU Lesser General Public License.  This exception does not
   however invalidate any other reasons why the executable file
   might be covered by the GNU Lesser General Public License.
   This exception applies to code released by its copyright holders
   in files containing the exception.  */


/* This file is needed by libio to define various configuration parameters.
   These are always the same in the GNU C library.  */


/* Define types for libio in terms of the standard internal type names.  */

/* stddef.h COPYRIGHT FUJITSU LIMITED 1998-2009                          */


 

 
 
 




/* Copyright (C) 1995-2004,2005,2006,2007 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */

/*
 *      ISO C99 Standard: 7.24
 *	Extended multibyte and wide character utilities	<wchar.h>
 */



/* stddef.h COPYRIGHT FUJITSU LIMITED 1998-2009                          */


 

 
 
 

typedef unsigned int wint_t;




/* We try to get wint_t from <stddef.h>, but not all GCC versions define it
   there.  So define it ourselves if it remains undefined.  */
/* Work around problems with the <stddef.h> file which doesn't put
   wint_t in the std namespace.  */

/* Conversion state information.  */
typedef struct
{
  int __count;
  union
  {
    wint_t __wch;
    char __wchb[4];
  } __value;		/* Value so far.  */
} __mbstate_t;


/* The rest of the file is only used if used if __need_mbstate_t is not
   defined.  */

typedef struct
{
  __off_t __pos;
  __mbstate_t __state;
} _G_fpos_t;
typedef struct
{
  __off64_t __pos;
  __mbstate_t __state;
} _G_fpos64_t;

typedef int _G_int16_t ;
typedef int _G_int32_t ;
typedef unsigned int _G_uint16_t ;
typedef unsigned int _G_uint32_t ;



/* These library features are always available in the GNU C library.  */




/* This is defined by <bits/stat.h> if `st_blksize' exists.  */


/* These are the vtbl details for ELF.  */



/* ALL of these should be defined in _G_config.h */

/* This define avoids name pollution if we're using GNU stdarg.h */
/* stdarg.h COPYRIGHT FUJITSU LIMITED 1998-2000                          */



typedef void *va_list;








/* For backward compatibility */




/* Magic numbers and bits for the _flags field.
   The magic numbers use the high-order bits of _flags;
   the remaining bits are available for variable flags.
   Note: The magic numbers must all be negative if stdio
   emulation is desired. */



/* These are "formatting flags" matching the iostream fmtflags enum values. */


struct _IO_jump_t;  struct _IO_FILE;

/* Handle lock.  */
typedef void _IO_lock_t;


/* A streammarker remembers a position in a buffer. */

struct _IO_marker {
  struct _IO_marker *_next;
  struct _IO_FILE *_sbuf;
  /* If _pos >= 0
 it points to _buf->Gbase()+_pos. FIXME comment */
  /* if _pos < 0, it points to _buf->eBptr()+_pos. FIXME comment */
  int _pos;
};

/* This is the structure from the libstdc++ codecvt class.  */
enum __codecvt_result
{
  __codecvt_ok,
  __codecvt_partial,
  __codecvt_error,
  __codecvt_noconv
};


struct _IO_FILE {
  int _flags;		/* High-order word is _IO_MAGIC; rest is flags. */

  /* The following pointers correspond to the C++ streambuf protocol. */
  /* Note:  Tk uses the _IO_read_ptr and _IO_read_end fields directly. */
  char* _IO_read_ptr;	/* Current read pointer */
  char* _IO_read_end;	/* End of get area. */
  char* _IO_read_base;	/* Start of putback+get area. */
  char* _IO_write_base;	/* Start of put area. */
  char* _IO_write_ptr;	/* Current put pointer. */
  char* _IO_write_end;	/* End of put area. */
  char* _IO_buf_base;	/* Start of reserve area. */
  char* _IO_buf_end;	/* End of reserve area. */
  /* The following fields are used to support backing up and undo. */
  char *_IO_save_base; /* Pointer to start of non-current get area. */
  char *_IO_backup_base;  /* Pointer to first valid character of backup area */
  char *_IO_save_end; /* Pointer to end of non-current get area. */

  struct _IO_marker *_markers;

  struct _IO_FILE *_chain;

  int _fileno;
  int _flags2;
  __off_t _old_offset; /* This used to be _offset but it's too small.  */

  /* 1+column number of pbase(); 0 is unknown. */
  unsigned short _cur_column;
  signed char _vtable_offset;
  char _shortbuf[1];

  /*  char* _save_gptr;  char* _save_egptr; */

  _IO_lock_t *_lock;
  __off64_t _offset;
  void *__pad1;
  void *__pad2;
  void *__pad3;
  void *__pad4;
  size_t __pad5;
  int _mode;
  /* Make sure we don't get into trouble again.  */
  char _unused2[15 * sizeof (int) - 4 * sizeof (void *) - sizeof (size_t)];
};

typedef struct _IO_FILE _IO_FILE;

struct _IO_FILE_plus;

extern struct _IO_FILE_plus _IO_2_1_stdin_;
extern struct _IO_FILE_plus _IO_2_1_stdout_;
extern struct _IO_FILE_plus _IO_2_1_stderr_;


/* Functions to do I/O and file management for a stream.  */

/* Read NBYTES bytes from COOKIE into a buffer pointed to by BUF.
   Return number of bytes read.  */
typedef __ssize_t __io_read_fn (void *__cookie, char *__buf, size_t __nbytes);

/* Write N bytes pointed to by BUF to COOKIE.  Write all N bytes
   unless there is an error.  Return number of bytes written, or -1 if
   there is an error without writing anything.  If the file has been
   opened for append (__mode.__append set), then set the file pointer
   to the end of the file and then do the write; if not, just write at
   the current file pointer.  */
typedef __ssize_t __io_write_fn (void *__cookie, const char *__buf,
				 size_t __n);

/* Move COOKIE's file position to *POS bytes from the
   beginning of the file (if W is SEEK_SET),
   the current position (if W is SEEK_CUR),
   or the end of the file (if W is SEEK_END).
   Set *POS to the new file position.
   Returns zero if successful, nonzero if not.  */
typedef int __io_seek_fn (void *__cookie, __off64_t *__pos, int __w);

/* Close COOKIE.  */
typedef int __io_close_fn (void *__cookie);





extern int __underflow (_IO_FILE *);
extern int __uflow (_IO_FILE *);
extern int __overflow (_IO_FILE *, int);





extern int _IO_getc (_IO_FILE *__fp);
extern int _IO_putc (int __c, _IO_FILE *__fp);
extern int _IO_feof (_IO_FILE *__fp) ;
extern int _IO_ferror (_IO_FILE *__fp) ;

extern int _IO_peekc_locked (_IO_FILE *__fp);

/* This one is for Emacs. */

extern void _IO_flockfile (_IO_FILE *) ;
extern void _IO_funlockfile (_IO_FILE *) ;
extern int _IO_ftrylockfile (_IO_FILE *) ;


extern int _IO_vfscanf (_IO_FILE * , const char * ,
			va_list, int *);
extern int _IO_vfprintf (_IO_FILE *, const char *,
			 va_list);
extern __ssize_t _IO_padn (_IO_FILE *, int, __ssize_t);
extern size_t _IO_sgetn (_IO_FILE *, void *, size_t);

extern __off64_t _IO_seekoff (_IO_FILE *, __off64_t, int, int);
extern __off64_t _IO_seekpos (_IO_FILE *, __off64_t, int);

extern void _IO_free_backup_area (_IO_FILE *) ;






/* The type of the second argument to `fgetpos' and `fsetpos'.  */

typedef _G_fpos_t fpos_t;


/* The possibilities for the third argument to `setvbuf'.  */


/* Default buffer size.  */


/* End of file character.
   Some things throughout the library rely on this being -1.  */


/* The possibilities for the third argument to `fseek'.
   These values should not be changed.  */


/* Default path prefix for `tempnam' and `tmpnam'.  */


/* Get the values:
   L_tmpnam	How long an array of chars must be to be passed to `tmpnam'.
   TMP_MAX	The minimum number of unique filenames generated by tmpnam
   		(and tempnam when it uses tmpnam's name space),
		or tempnam (the two are separate).
   L_ctermid	How long an array to pass to `ctermid'.
   L_cuserid	How long an array to pass to `cuserid'.
   FOPEN_MAX	Minimum number of files that can be open at once.
   FILENAME_MAX	Maximum length of a filename.  */
/* Copyright (C) 1994, 1997, 1998, 1999 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */







/* Standard streams.  */
extern struct _IO_FILE *stdin;		/* Standard input stream.  */
extern struct _IO_FILE *stdout;		/* Standard output stream.  */
extern struct _IO_FILE *stderr;		/* Standard error output stream.  */
/* C89/C99 say they're macros.  Make them happy.  */


/* Remove file FILENAME.  */
extern int remove (const char *__filename) ;
/* Rename file OLD to NEW.  */
extern int rename (const char *__old, const char *__new) ;




/* Create a temporary file and open it read/write.

   This function is a possible cancellation points and therefore not
   marked with __THROW.  */
extern FILE *tmpfile (void) ;


/* Generate a temporary filename.  */
extern char *tmpnam (char *__s)  ;


/* This is the reentrant variant of `tmpnam'.  The only difference is
   that it does not allow S to be NULL.  */
extern char *tmpnam_r (char *__s)  ;


/* Generate a unique temporary filename using up to five characters of PFX
   if it is not NULL.  The directory to put this file in is searched for
   as follows: First the environment variable "TMPDIR" is checked.
   If it contains the name of a writable directory, that directory is used.
   If not and if DIR is not NULL, that value is checked.  If that fails,
   P_tmpdir is tried and finally "/tmp".  The storage for the filename
   is allocated by `malloc'.  */
extern char *tempnam (const char *__dir, const char *__pfx)
       ;



/* Close STREAM.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern int fclose (FILE *__stream);
/* Flush STREAM, or all streams if STREAM is NULL.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern int fflush (FILE *__stream);


/* Faster versions when locking is not required.

   This function is not part of POSIX and therefore no official
   cancellation point.  But due to similarity with an POSIX interface
   or due to the implementation it is a cancellation point and
   therefore not marked with __THROW.  */
extern int fflush_unlocked (FILE *__stream);




/* Open a file and create a new stream for it.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern FILE *fopen (const char * __filename,
		    const char * __modes) ;
/* Open a file, replacing an existing stream with it.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern FILE *freopen (const char * __filename,
		      const char * __modes,
		      FILE * __stream) ;


/* Create a new stream that refers to an existing system file descriptor.  */
extern FILE *fdopen (int __fd, const char *__modes)  ;




/* If BUF is NULL, make STREAM unbuffered.
   Else make it use buffer BUF, of size BUFSIZ.  */
extern void setbuf (FILE * __stream, char * __buf) ;
/* Make STREAM use buffering mode MODE.
   If BUF is not NULL, use N bytes of it for buffering;
   else allocate an internal buffer N bytes long.  */
extern int setvbuf (FILE * __stream, char * __buf,
		    int __modes, size_t __n) ;


/* If BUF is NULL, make STREAM unbuffered.
   Else make it use SIZE bytes of BUF for buffering.  */
extern void setbuffer (FILE * __stream, char * __buf,
		       size_t __size) ;

/* Make STREAM line-buffered.  */
extern void setlinebuf (FILE *__stream) ;



/* Write formatted output to STREAM.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern int fprintf (FILE * __stream,
		    const char * __format, ...);
/* Write formatted output to stdout.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern int printf (const char * __format, ...);
/* Write formatted output to S.  */
extern int sprintf (char * __s,
		    const char * __format, ...) ;

/* Write formatted output to S from argument list ARG.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern int vfprintf (FILE * __s, const char * __format,
		     va_list __arg);
/* Write formatted output to stdout from argument list ARG.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern int vprintf (const char * __format, va_list __arg);
/* Write formatted output to S from argument list ARG.  */
extern int vsprintf (char * __s, const char * __format,
		     va_list __arg) ;



/* Maximum chars of output to write in MAXLEN.  */
extern int snprintf (char * __s, size_t __maxlen,
		     const char * __format, ...)
      ;

extern int vsnprintf (char * __s, size_t __maxlen,
		      const char * __format, va_list __arg)
      ;





/* Read formatted input from STREAM.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern int fscanf (FILE * __stream,
		   const char * __format, ...) ;
/* Read formatted input from stdin.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern int scanf (const char * __format, ...) ;
/* Read formatted input from S.  */
extern int sscanf (const char * __s,
		   const char * __format, ...) ;

extern int __isoc99_fscanf (FILE * __stream,
			    const char * __format, ...) ;
extern int __isoc99_scanf (const char * __format, ...) ;
extern int __isoc99_sscanf (const char * __s,
			    const char * __format, ...) ;




/* Read formatted input from S into argument list ARG.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern int vfscanf (FILE * __s, const char * __format,
		    va_list __arg)
      ;

/* Read formatted input from stdin into argument list ARG.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern int vscanf (const char * __format, va_list __arg)
      ;

/* Read formatted input from S into argument list ARG.  */
extern int vsscanf (const char * __s,
		    const char * __format, va_list __arg)
      ;

extern int __isoc99_vfscanf (FILE * __s,
			     const char * __format,
			     va_list __arg) ;
extern int __isoc99_vscanf (const char * __format,
			    va_list __arg) ;
extern int __isoc99_vsscanf (const char * __s,
			     const char * __format,
			     va_list __arg) ;





/* Read a character from STREAM.

   These functions are possible cancellation points and therefore not
   marked with __THROW.  */
extern int fgetc (FILE *__stream);
extern int getc (FILE *__stream);

/* Read a character from stdin.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern int getchar (void);


/* The C standard explicitly says this is a macro, so we always do the
   optimization for it.  */

/* These are defined in POSIX.1:1996.

   These functions are possible cancellation points and therefore not
   marked with __THROW.  */
extern int getc_unlocked (FILE *__stream);
extern int getchar_unlocked (void);

/* Faster version when locking is not necessary.

   This function is not part of POSIX and therefore no official
   cancellation point.  But due to similarity with an POSIX interface
   or due to the implementation it is a cancellation point and
   therefore not marked with __THROW.  */
extern int fgetc_unlocked (FILE *__stream);



/* Write a character to STREAM.

   These functions are possible cancellation points and therefore not
   marked with __THROW.

   These functions is a possible cancellation point and therefore not
   marked with __THROW.  */
extern int fputc (int __c, FILE *__stream);
extern int putc (int __c, FILE *__stream);

/* Write a character to stdout.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern int putchar (int __c);


/* The C standard explicitly says this can be a macro,
   so we always do the optimization for it.  */

/* Faster version when locking is not necessary.

   This function is not part of POSIX and therefore no official
   cancellation point.  But due to similarity with an POSIX interface
   or due to the implementation it is a cancellation point and
   therefore not marked with __THROW.  */
extern int fputc_unlocked (int __c, FILE *__stream);

/* These are defined in POSIX.1:1996.

   These functions are possible cancellation points and therefore not
   marked with __THROW.  */
extern int putc_unlocked (int __c, FILE *__stream);
extern int putchar_unlocked (int __c);


/* Get a word (int) from STREAM.  */
extern int getw (FILE *__stream);

/* Write a word (int) to STREAM.  */
extern int putw (int __w, FILE *__stream);



/* Get a newline-terminated string of finite length from STREAM.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern char *fgets (char * __s, int __n, FILE * __stream)
     ;

/* Get a newline-terminated string from stdin, removing the newline.
   DO NOT USE THIS FUNCTION!!  There is no limit on how much it will read.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern char *gets (char *__s) ;







/* Write a string to STREAM.

   This function is a possible cancellation points and therefore not
   marked with __THROW.  */
extern int fputs (const char * __s, FILE * __stream);

/* Write a string, followed by a newline, to stdout.

   This function is a possible cancellation points and therefore not
   marked with __THROW.  */
extern int puts (const char *__s);


/* Push a character back onto the input buffer of STREAM.

   This function is a possible cancellation points and therefore not
   marked with __THROW.  */
extern int ungetc (int __c, FILE *__stream);


/* Read chunks of generic data from STREAM.

   This function is a possible cancellation points and therefore not
   marked with __THROW.  */
extern size_t fread (void * __ptr, size_t __size,
		     size_t __n, FILE * __stream) ;
/* Write chunks of generic data to STREAM.

   This function is a possible cancellation points and therefore not
   marked with __THROW.  */
extern size_t fwrite (const void * __ptr, size_t __size,
		      size_t __n, FILE * __s) ;



/* Faster versions when locking is not necessary.

   These functions are not part of POSIX and therefore no official
   cancellation point.  But due to similarity with an POSIX interface
   or due to the implementation they are cancellation points and
   therefore not marked with __THROW.  */
extern size_t fread_unlocked (void * __ptr, size_t __size,
			      size_t __n, FILE * __stream) ;
extern size_t fwrite_unlocked (const void * __ptr, size_t __size,
			       size_t __n, FILE * __stream) ;



/* Seek to a certain position on STREAM.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern int fseek (FILE *__stream, long int __off, int __whence);
/* Return the current position of STREAM.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern long int ftell (FILE *__stream) ;
/* Rewind to the beginning of STREAM.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern void rewind (FILE *__stream);


/* The Single Unix Specification, Version 2, specifies an alternative,
   more adequate interface for the two functions above which deal with
   file offset.  `long int' is not the right type.  These definitions
   are originally defined in the Large File Support API.  */

/* Seek to a certain position on STREAM.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern int fseeko (FILE *__stream, __off_t __off, int __whence);
/* Return the current position of STREAM.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern __off_t ftello (FILE *__stream) ;


/* Get STREAM's position.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern int fgetpos (FILE * __stream, fpos_t * __pos);
/* Set STREAM's position.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern int fsetpos (FILE *__stream, const fpos_t *__pos);




/* Clear the error and EOF indicators for STREAM.  */
extern void clearerr (FILE *__stream) ;
/* Return the EOF indicator for STREAM.  */
extern int feof (FILE *__stream)  ;
/* Return the error indicator for STREAM.  */
extern int ferror (FILE *__stream)  ;


/* Faster versions when locking is not required.  */
extern void clearerr_unlocked (FILE *__stream) ;
extern int feof_unlocked (FILE *__stream)  ;
extern int ferror_unlocked (FILE *__stream)  ;



/* Print a message describing the meaning of the value of errno.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern void perror (const char *__s);


/* Provide the declarations for `sys_errlist' and `sys_nerr' if they
   are available on this system.  Even if available, these variables
   should not be used directly.  The `strerror' function provides
   all the necessary functionality.  */
/* Declare sys_errlist and sys_nerr, or don't.  Compatibility (do) version.
   Copyright (C) 2002 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */


/* sys_errlist and sys_nerr are deprecated.  Use strerror instead.  */

extern int sys_nerr;
extern const char *const sys_errlist[];


/* Return the system file descriptor for STREAM.  */
extern int fileno (FILE *__stream)  ;

/* Faster version when locking is not required.  */
extern int fileno_unlocked (FILE *__stream)  ;


/* Create a new stream connected to a pipe running the given command.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern FILE *popen (const char *__command, const char *__modes) ;

/* Close a stream opened by popen and return the status of its child.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern int pclose (FILE *__stream);


/* Return the name of the controlling terminal.  */
extern char *ctermid (char *__s) ;






/* These are defined in POSIX.1:1996.  */

/* Acquire ownership of STREAM.  */
extern void flockfile (FILE *__stream) ;

/* Try to acquire ownership of STREAM but do not block if it is not
   possible.  */
extern int ftrylockfile (FILE *__stream)  ;

/* Relinquish the ownership granted for STREAM.  */
extern void funlockfile (FILE *__stream) ;


/* If we are compiling with optimizing read this file.  It contains
   several optimizing inline functions and macros.  */




/* Copyright (C) 1991-2003,2004,2005,2006,2007 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */

/*
 *	ISO C99 Standard: 7.20 General utilities	<stdlib.h>
 */


/* features.h Copyright FUJITSU LIMITED 2008-2010 */


/* Get size_t, wchar_t and NULL from <stddef.h>.  */
/* stddef.h COPYRIGHT FUJITSU LIMITED 1998-2009                          */


 

 
 
 
typedef int wchar_t;










/* Returned by `div'.  */
typedef struct
  {
    int quot;			/* Quotient.  */
    int rem;			/* Remainder.  */
  } div_t;

/* Returned by `ldiv'.  */
typedef struct
  {
    long int quot;		/* Quotient.  */
    long int rem;		/* Remainder.  */
  } ldiv_t;



/* Returned by `lldiv'.  */
 typedef struct
  {
    long long int quot;		/* Quotient.  */
    long long int rem;		/* Remainder.  */
  } lldiv_t;



/* The largest number rand will return (same as INT_MAX).  */


/* We define these the same for all machines.
   Changes from this to the outside world should be done in `_exit'.  */


/* Maximum length of a multibyte character in the current locale.  */
extern size_t __ctype_get_mb_cur_max (void)  ;



/* Convert a string to a floating-point number.  */
extern double atof (const char *__nptr)
        ;
/* Convert a string to an integer.  */
extern int atoi (const char *__nptr)
        ;
/* Convert a string to a long integer.  */
extern long int atol (const char *__nptr)
        ;



/* Convert a string to a long long integer.  */
 extern long long int atoll (const char *__nptr)
        ;



/* Convert a string to a floating-point number.  */
extern double strtod (const char * __nptr,
		      char ** __endptr)
       ;



/* Likewise for `float' and `long double' sizes of floating-point numbers.  */
extern float strtof (const char * __nptr,
		     char ** __endptr)   ;

extern long double strtold (const char * __nptr,
			    char ** __endptr)
       ;



/* Convert a string to a long integer.  */
extern long int strtol (const char * __nptr,
			char ** __endptr, int __base)
       ;
/* Convert a string to an unsigned long integer.  */
extern unsigned long int strtoul (const char * __nptr,
				  char ** __endptr, int __base)
       ;


/* Convert a string to a quadword integer.  */

extern long long int strtoq (const char * __nptr,
			     char ** __endptr, int __base)
       ;
/* Convert a string to an unsigned quadword integer.  */

extern unsigned long long int strtouq (const char * __nptr,
				       char ** __endptr, int __base)
       ;


/* Convert a string to a quadword integer.  */

extern long long int strtoll (const char * __nptr,
			      char ** __endptr, int __base)
       ;
/* Convert a string to an unsigned quadword integer.  */

extern unsigned long long int strtoull (const char * __nptr,
					char ** __endptr, int __base)
       ;







/* Convert N to base 64 using the digits "./0-9A-Za-z", least-significant
   digit first.  Returns a pointer to static storage overwritten by the
   next call.  */
extern char *l64a (long int __n)  ;

/* Read a number from a string S in base 64 as above.  */
extern long int a64l (const char *__s)
        ;


/* Copyright (C) 1991,1992,1994,1995,1996,1997,1998,1999,2000,2001,2002,2006
   	Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */

/*
 *	POSIX Standard: 2.6 Primitive System Data Types	<sys/types.h>
 */


/* features.h Copyright FUJITSU LIMITED 2008-2010 */





typedef __u_char u_char;
typedef __u_short u_short;
typedef __u_int u_int;
typedef __u_long u_long;
typedef __quad_t quad_t;
typedef __u_quad_t u_quad_t;
typedef __fsid_t fsid_t;

typedef __loff_t loff_t;

typedef __ino_t ino_t;

typedef __dev_t dev_t;

typedef __gid_t gid_t;

typedef __mode_t mode_t;

typedef __nlink_t nlink_t;

typedef __uid_t uid_t;

typedef __off_t off_t;

typedef __pid_t pid_t;

typedef __id_t id_t;

typedef __ssize_t ssize_t;

typedef __daddr_t daddr_t;
typedef __caddr_t caddr_t;

typedef __key_t key_t;

/* Copyright (C) 1991-1999,2000,2001,2002,2003,2006
	Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */

/*
 *	ISO C99 Standard: 7.23 Date and time	<time.h>
 */








/* Returned by `time'.  */
typedef __time_t time_t;






/* Clock ID used in clock and timer functions.  */
typedef __clockid_t clockid_t;




/* Timer ID returned by `timer_create'.  */
typedef __timer_t timer_t;








/* stddef.h COPYRIGHT FUJITSU LIMITED 1998-2009                          */


 

 
 
 





/* Old compatibility names for C types.  */
typedef unsigned long int ulong;
typedef unsigned short int ushort;
typedef unsigned int uint;

/* These size-specific names are used by some of the inet code.  */


/* These types are defined by the ISO C99 header <inttypes.h>. */
typedef	char int8_t;
typedef	short int int16_t;
typedef	int int32_t;
typedef long int int64_t;

/* But these were defined by ISO C without the first `_'.  */
typedef	unsigned char u_int8_t;
typedef	unsigned short int u_int16_t;
typedef	unsigned int u_int32_t;
typedef unsigned long int u_int64_t;

typedef int register_t;



/* In BSD <sys/types.h> is expected to define BYTE_ORDER.  */
/* Copyright (C) 1992, 1996, 1997, 2000 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */


/* features.h Copyright FUJITSU LIMITED 2008-2010 */


/* Definitions for byte order, according to significance of bytes,
   from low addresses to high addresses.  The value is what you get by
   putting '4' in the most significant byte, '3' in the second most
   significant byte, '2' in the second least significant byte, and '1'
   in the least significant byte, and then writing down one digit for
   each byte, starting with the byte at the lowest address at the left,
   and proceeding to the byte with the highest address at the right.  */


/* This file defines `__BYTE_ORDER' for the particular machine.  */
/* Sparc is big-endian, but v9 supports endian conversion on loads/stores
   and GCC supports such a mode.  Be prepared.  */



/* Some machines may need to use a different endianness for floating point
   values.  */




/* It also defines `fd_set' and the FD_* macros for `select'.  */
/* `fd_set' type and related macros, and `select'/`pselect' declarations.
   Copyright (C) 1996,97,98,99,2000,01,02,2003 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */

/*	POSIX 1003.1g: 6.2 Select from File Descriptor Sets <sys/select.h>  */


/* features.h Copyright FUJITSU LIMITED 2008-2010 */


/* Get definition of needed basic types.  */

/* Get __FD_* definitions.  */
/* Copyright (C) 1997, 1998, 2001 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */



/* We don't use `memset' because this would require a prototype and
   the array isn't too big.  */

/* Get __sigset_t.  */
/* __sig_atomic_t, __sigset_t, and related definitions.  Linux version.
   Copyright (C) 1991, 1992, 1994, 1996, 1997, 2007
   Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */


typedef int __sig_atomic_t;

/* A `sigset_t' has a bit for each signal.  */

typedef struct
  {
    unsigned long int __val[(1024 / (8 * sizeof (unsigned long int)))];
  } __sigset_t;



/* We only want to define these functions if <signal.h> was actually
   included; otherwise we were included just to define the types.  Since we
   are namespace-clean, it wouldn't hurt to define extra macros.  But
   trouble can be caused by functions being defined (e.g., any global
   register vars declared later will cause compilation errors).  */


typedef __sigset_t sigset_t;

/* Get definition of timer specification structures.  */
/* Copyright (C) 1991-1999,2000,2001,2002,2003,2006
	Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */

/*
 *	ISO C99 Standard: 7.23 Date and time	<time.h>
 */











/* POSIX.1b structure for a time value.  This is like a `struct timeval' but
   has nanoseconds instead of microseconds.  */
struct timespec
  {
    __time_t tv_sec;		/* Seconds.  */
    long int tv_nsec;		/* Nanoseconds.  */
  };




/* System-dependent timing definitions.  Generic version.
   Copyright (C) 1996,1997,1999-2002,2003 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */

/*
 * Never include this file directly; use <time.h> instead.
 */



/* A time value that is accurate to the nearest
   microsecond but also has a range of years.  */
struct timeval
  {
    __time_t tv_sec;		/* Seconds.  */
    __suseconds_t tv_usec;	/* Microseconds.  */
  };

typedef __suseconds_t suseconds_t;


/* The fd_set member is required to be an array of longs.  */
typedef long int __fd_mask;

/* Some versions of <linux/posix_types.h> define these macros.  */
/* It's easier to assume 8-bit bytes than to get CHAR_BIT.  */

/* fd_set for select and pselect.  */
typedef struct
  {
    /* XPG4.2 requires this member name.  Otherwise avoid the name
       from the global namespace.  */
    __fd_mask __fds_bits[1024 / (8 * sizeof (__fd_mask))];
  } fd_set;

/* Maximum number of file descriptors in `fd_set'.  */

/* Sometimes the fd_set member is assumed to have this type.  */
typedef __fd_mask fd_mask;

/* Number of bits per word of `fd_set' (some code assumes this is 32).  */


/* Access macros for `fd_set'.  */




/* Check the first NFDS descriptors each in READFDS (if not NULL) for read
   readiness, in WRITEFDS (if not NULL) for write readiness, and in EXCEPTFDS
   (if not NULL) for exceptional conditions.  If TIMEOUT is not NULL, time out
   after waiting the interval specified therein.  Returns the number of ready
   descriptors, or -1 for errors.

   This function is a cancellation point and therefore not marked with
   __THROW.  */
extern int select (int __nfds, fd_set * __readfds,
		   fd_set * __writefds,
		   fd_set * __exceptfds,
		   struct timeval * __timeout);

/* Same as above only that the TIMEOUT value is given with higher
   resolution and a sigmask which is been set temporarily.  This version
   should be used.

   This function is a cancellation point and therefore not marked with
   __THROW.  */
extern int pselect (int __nfds, fd_set * __readfds,
		    fd_set * __writefds,
		    fd_set * __exceptfds,
		    const struct timespec * __timeout,
		    const __sigset_t * __sigmask);




/* BSD defines these symbols, so we follow.  */
/* Definitions of macros to access `dev_t' values.
   Copyright (C) 1996, 1997, 1999, 2003, 2004, 2007
   Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */


/* features.h Copyright FUJITSU LIMITED 2008-2010 */


/* If the compiler does not know long long it is out of luck.  We are
   not going to hack weird hacks to support the dev_t representation
   they need.  */

extern unsigned int gnu_dev_major (unsigned long long int __dev)
     ;

extern unsigned int gnu_dev_minor (unsigned long long int __dev)
     ;

extern unsigned long long int gnu_dev_makedev (unsigned int __major,
					       unsigned int __minor)
     ;



/* Access the functions with their traditional names.  */




/* Types from the Large File Support interface.  */
typedef __blkcnt_t blkcnt_t;	 /* Type to count number of disk blocks.  */
typedef __fsblkcnt_t fsblkcnt_t; /* Type to count file system blocks.  */
typedef __fsfilcnt_t fsfilcnt_t; /* Type to count file system inodes.  */



/* Now add the thread types.  */
/* Machine-specific pthread type layouts.  SPARC version.
   Copyright (C) 2003, 2004, 2005, 2006, 2007 Free Software Foundation, Inc.
   This file is part of the GNU C Library.
   Contributed by Jakub Jelinek <jakub@redhat.com>, 2003.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */


/* Determine the wordsize from the preprocessor defines.  */






/* Thread identifiers.  The structure of the attribute type is
   deliberately not exposed.  */
typedef unsigned long int pthread_t;


typedef union
{
  char __size[56];
  long int __align;
} pthread_attr_t;


typedef struct __pthread_internal_list
{
  struct __pthread_internal_list *__prev;
  struct __pthread_internal_list *__next;
} __pthread_list_t;


/* Data structures for mutex handling.  The structure of the attribute
   type is deliberately not exposed.  */
typedef union
{
  struct __pthread_mutex_s
  {
    int __lock;
    unsigned int __count;
    int __owner;
    unsigned int __nusers;
    /* KIND must stay at this position in the structure to maintain
       binary compatibility.  */
    int __kind;
    int __spins;
    __pthread_list_t __list;
  } __data;
  char __size[40];
  long int __align;
} pthread_mutex_t;

typedef union
{
  char __size[4];
  int __align;
} pthread_mutexattr_t;


/* Data structure for conditional variable handling.  The structure of
   the attribute type is deliberately not exposed.  */
typedef union
{
  struct
  {
    int __lock;
    unsigned int __futex;
     unsigned long long int __total_seq;
     unsigned long long int __wakeup_seq;
     unsigned long long int __woken_seq;
    void *__mutex;
    unsigned int __nwaiters;
    unsigned int __broadcast_seq;
  } __data;
  char __size[48];
   long long int __align;
} pthread_cond_t;

typedef union
{
  char __size[4];
  int __align;
} pthread_condattr_t;


/* Keys for thread-specific data */
typedef unsigned int pthread_key_t;


/* Once-only execution */
typedef int pthread_once_t;


/* Data structure for read-write lock variable handling.  The
   structure of the attribute type is deliberately not exposed.  */
typedef union
{
  struct
  {
    int __lock;
    unsigned int __nr_readers;
    unsigned int __readers_wakeup;
    unsigned int __writer_wakeup;
    unsigned int __nr_readers_queued;
    unsigned int __nr_writers_queued;
    int __writer;
    int __shared;
    unsigned long int __pad1;
    unsigned long int __pad2;
    /* FLAGS must stay at this position in the structure to maintain
       binary compatibility.  */
    unsigned int __flags;
  } __data;
  char __size[56];
  long int __align;
} pthread_rwlock_t;

typedef union
{
  char __size[8];
  long int __align;
} pthread_rwlockattr_t;


/* POSIX spinlock data type.  */
typedef volatile int pthread_spinlock_t;


/* POSIX barriers data type.  The structure of the type is
   deliberately not exposed.  */
typedef union
{
  char __size[32];
  long int __align;
} pthread_barrier_t;

typedef union
{
  char __size[4];
  int __align;
} pthread_barrierattr_t;






/* These are the functions that actually do things.  The `random', `srandom',
   `initstate' and `setstate' functions are those from BSD Unices.
   The `rand' and `srand' functions are required by the ANSI standard.
   We provide both interfaces to the same random number generator.  */
/* Return a random long integer between 0 and RAND_MAX inclusive.  */
extern long int random (void) ;

/* Seed the random number generator with the given number.  */
extern void srandom (unsigned int __seed) ;

/* Initialize the random number generator to use state buffer STATEBUF,
   of length STATELEN, and seed it with SEED.  Optimal lengths are 8, 16,
   32, 64, 128 and 256, the bigger the better; values less than 8 will
   cause an error and values greater than 256 will be rounded down.  */
extern char *initstate (unsigned int __seed, char *__statebuf,
			size_t __statelen)  ;

/* Switch the random number generator to state buffer STATEBUF,
   which should have been previously initialized by `initstate'.  */
extern char *setstate (char *__statebuf)  ;


/* Reentrant versions of the `random' family of functions.
   These functions all use the following data structure to contain
   state, rather than global state variables.  */

struct random_data
  {
    int32_t *fptr;		/* Front pointer.  */
    int32_t *rptr;		/* Rear pointer.  */
    int32_t *state;		/* Array of state values.  */
    int rand_type;		/* Type of random number generator.  */
    int rand_deg;		/* Degree of random number generator.  */
    int rand_sep;		/* Distance between front and rear.  */
    int32_t *end_ptr;		/* Pointer behind state table.  */
  };

extern int random_r (struct random_data * __buf,
		     int32_t * __result)  ;

extern int srandom_r (unsigned int __seed, struct random_data *__buf)
      ;

extern int initstate_r (unsigned int __seed, char * __statebuf,
			size_t __statelen,
			struct random_data * __buf)
      ;

extern int setstate_r (char * __statebuf,
		       struct random_data * __buf)
      ;



/* Return a random integer between 0 and RAND_MAX inclusive.  */
extern int rand (void) ;
/* Seed the random number generator with the given number.  */
extern void srand (unsigned int __seed) ;


/* Reentrant interface according to POSIX.1.  */
extern int rand_r (unsigned int *__seed) ;


/* System V style 48-bit random number generator functions.  */

/* Return non-negative, double-precision floating-point value in [0.0,1.0).  */
extern double drand48 (void) ;
extern double erand48 (unsigned short int __xsubi[3])  ;

/* Return non-negative, long integer in [0,2^31).  */
extern long int lrand48 (void) ;
extern long int nrand48 (unsigned short int __xsubi[3])
      ;

/* Return signed, long integers in [-2^31,2^31).  */
extern long int mrand48 (void) ;
extern long int jrand48 (unsigned short int __xsubi[3])
      ;

/* Seed random number generator.  */
extern void srand48 (long int __seedval) ;
extern unsigned short int *seed48 (unsigned short int __seed16v[3])
      ;
extern void lcong48 (unsigned short int __param[7])  ;

/* Data structure for communication with thread safe versions.  This
   type is to be regarded as opaque.  It's only exported because users
   have to allocate objects of this type.  */
struct drand48_data
  {
    unsigned short int __x[3];	/* Current state.  */
    unsigned short int __old_x[3]; /* Old state.  */
    unsigned short int __c;	/* Additive const. in congruential formula.  */
    unsigned short int __init;	/* Flag for initializing.  */
    unsigned long long int __a;	/* Factor in congruential formula.  */
  };

/* Return non-negative, double-precision floating-point value in [0.0,1.0).  */
extern int drand48_r (struct drand48_data * __buffer,
		      double * __result)  ;
extern int erand48_r (unsigned short int __xsubi[3],
		      struct drand48_data * __buffer,
		      double * __result)  ;

/* Return non-negative, long integer in [0,2^31).  */
extern int lrand48_r (struct drand48_data * __buffer,
		      long int * __result)
      ;
extern int nrand48_r (unsigned short int __xsubi[3],
		      struct drand48_data * __buffer,
		      long int * __result)
      ;

/* Return signed, long integers in [-2^31,2^31).  */
extern int mrand48_r (struct drand48_data * __buffer,
		      long int * __result)
      ;
extern int jrand48_r (unsigned short int __xsubi[3],
		      struct drand48_data * __buffer,
		      long int * __result)
      ;

/* Seed random number generator.  */
extern int srand48_r (long int __seedval, struct drand48_data *__buffer)
      ;

extern int seed48_r (unsigned short int __seed16v[3],
		     struct drand48_data *__buffer)  ;

extern int lcong48_r (unsigned short int __param[7],
		      struct drand48_data *__buffer)
      ;



/* Allocate SIZE bytes of memory.  */
extern void *malloc (size_t __size)   ;
/* Allocate NMEMB elements of SIZE bytes each, all initialized to 0.  */
extern void *calloc (size_t __nmemb, size_t __size)
       ;



/* Re-allocate the previously allocated block
   in PTR, making the new block SIZE bytes long.  */
/* __attribute_malloc__ is not used, because if realloc returns
   the same pointer that was passed to it, aliasing needs to be allowed
   between objects pointed by the old and new pointers.  */
extern void *realloc (void *__ptr, size_t __size)
      ;
/* Free a block allocated by `malloc', `realloc' or `calloc'.  */
extern void free (void *__ptr) ;


/* Free a block.  An alias for `free'.	(Sun Unices).  */
extern void cfree (void *__ptr) ;

/* Copyright (C) 1992, 1996, 1997, 1998, 1999 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */


/* features.h Copyright FUJITSU LIMITED 2008-2010 */


/* stddef.h COPYRIGHT FUJITSU LIMITED 1998-2009                          */


 

 
 
 







/* Remove any previous definitions.  */

/* Allocate a block that will be freed when the calling function exits.  */
extern void *alloca (size_t __size) ;





/* Allocate SIZE bytes on a page boundary.  The storage cannot be freed.  */
extern void *valloc (size_t __size)   ;

/* Allocate memory of SIZE bytes with an alignment of ALIGNMENT.  */
extern int posix_memalign (void **__memptr, size_t __alignment, size_t __size)
       ;


/* Abort execution and generate a core-dump.  */
extern void abort (void)  ;


/* Register a function to be called when `exit' is called.  */
extern int atexit (void (*__func) (void))  ;


/* Register a function to be called with the status
   given to `exit' and the given argument.  */
extern int on_exit (void (*__func) (int __status, void *__arg), void *__arg)
      ;


/* Call all functions registered with `atexit' and `on_exit',
   in the reverse of the order in which they were registered
   perform stdio cleanup, and terminate program execution with STATUS.  */
extern void exit (int __status)  ;



/* Terminate the program with STATUS without calling any of the
   functions registered with `atexit' or `on_exit'.  */
extern void _Exit (int __status)  ;




/* Return the value of envariable NAME, or NULL if it doesn't exist.  */
extern char *getenv (const char *__name)   ;


/* This function is similar to the above but returns NULL if the
   programs is running with SUID or SGID enabled.  */
extern char *__secure_getenv (const char *__name)
       ;

/* The SVID says this is in <stdio.h>, but this seems a better place.	*/
/* Put STRING, which is of the form "NAME=VALUE", in the environment.
   If there is no `=', remove NAME from the environment.  */
extern int putenv (char *__string)  ;

/* Set NAME to VALUE in the environment.
   If REPLACE is nonzero, overwrite an existing value.  */
extern int setenv (const char *__name, const char *__value, int __replace)
      ;

/* Remove the variable NAME from the environment.  */
extern int unsetenv (const char *__name) ;

/* The `clearenv' was planned to be added to POSIX.1 but probably
   never made it.  Nevertheless the POSIX.9 standard (POSIX bindings
   for Fortran 77) requires this function.  */
extern int clearenv (void) ;


/* Generate a unique temporary file name from TEMPLATE.
   The last six characters of TEMPLATE must be "XXXXXX";
   they are replaced with a string that makes the file name unique.
   Returns TEMPLATE, or a null pointer if it cannot get a unique file name.  */
extern char *mktemp (char *__template)   ;

/* Generate a unique temporary file name from TEMPLATE.
   The last six characters of TEMPLATE must be "XXXXXX";
   they are replaced with a string that makes the filename unique.
   Returns a file descriptor open on the file for reading and writing,
   or -1 if it cannot create a uniquely-named file.

   This function is a possible cancellation points and therefore not
   marked with __THROW.  */
extern int mkstemp (char *__template)  ;

/* Create a unique temporary directory from TEMPLATE.
   The last six characters of TEMPLATE must be "XXXXXX";
   they are replaced with a string that makes the directory name unique.
   Returns TEMPLATE, or a null pointer if it cannot get a unique name.
   The directory is created mode 700.  */
extern char *mkdtemp (char *__template)   ;




/* Execute the given line as a shell command.

   This function is a cancellation point and therefore not marked with
   __THROW.  */
extern int system (const char *__command) ;




/* Return the canonical absolute name of file NAME.  If RESOLVED is
   null, the result is malloc'd; otherwise, if the canonical name is
   PATH_MAX chars or more, returns null with `errno' set to
   ENAMETOOLONG; if the name fits in fewer than PATH_MAX chars,
   returns the name in RESOLVED.  */
extern char *realpath (const char * __name,
		       char * __resolved)  ;


/* Shorthand for type of comparison functions.  */
typedef int (*__compar_fn_t) (const void *, const void *);



/* Do a binary search for KEY in BASE, which consists of NMEMB elements
   of SIZE bytes each, using COMPAR to perform the comparisons.  */
extern void *bsearch (const void *__key, const void *__base,
		      size_t __nmemb, size_t __size, __compar_fn_t __compar)
      ;

/* Sort NMEMB elements of BASE, of SIZE bytes each,
   using COMPAR to perform the comparisons.  */
extern void qsort (void *__base, size_t __nmemb, size_t __size,
		   __compar_fn_t __compar) ;


/* Return the absolute value of X.  */
extern int abs (int __x)   ;
extern long int labs (long int __x)   ;


 extern long long int llabs (long long int __x)
       ;



/* Return the `div_t', `ldiv_t' or `lldiv_t' representation
   of the value of NUMER over DENOM. */
/* GCC may have built-ins for these someday.  */
extern div_t div (int __numer, int __denom)
       ;
extern ldiv_t ldiv (long int __numer, long int __denom)
       ;



 extern lldiv_t lldiv (long long int __numer,
				    long long int __denom)
       ;



/* Convert floating point numbers to strings.  The returned values are
   valid only until another call to the same function.  */

/* Convert VALUE to a string with NDIGIT digits and return a pointer to
   this.  Set *DECPT with the position of the decimal character and *SIGN
   with the sign of the number.  */
extern char *ecvt (double __value, int __ndigit, int * __decpt,
		   int * __sign)   ;

/* Convert VALUE to a string rounded to NDIGIT decimal digits.  Set *DECPT
   with the position of the decimal character and *SIGN with the sign of
   the number.  */
extern char *fcvt (double __value, int __ndigit, int * __decpt,
		   int * __sign)   ;

/* If possible convert VALUE to a string with NDIGIT significant digits.
   Otherwise use exponential representation.  The resulting string will
   be written to BUF.  */
extern char *gcvt (double __value, int __ndigit, char *__buf)
       ;


/* Long double versions of above functions.  */
extern char *qecvt (long double __value, int __ndigit,
		    int * __decpt, int * __sign)
       ;
extern char *qfcvt (long double __value, int __ndigit,
		    int * __decpt, int * __sign)
       ;
extern char *qgcvt (long double __value, int __ndigit, char *__buf)
       ;


/* Reentrant version of the functions above which provide their own
   buffers.  */
extern int ecvt_r (double __value, int __ndigit, int * __decpt,
		   int * __sign, char * __buf,
		   size_t __len)  ;
extern int fcvt_r (double __value, int __ndigit, int * __decpt,
		   int * __sign, char * __buf,
		   size_t __len)  ;

extern int qecvt_r (long double __value, int __ndigit,
		    int * __decpt, int * __sign,
		    char * __buf, size_t __len)
      ;
extern int qfcvt_r (long double __value, int __ndigit,
		    int * __decpt, int * __sign,
		    char * __buf, size_t __len)
      ;



/* Return the length of the multibyte character
   in S, which is no longer than N.  */
extern int mblen (const char *__s, size_t __n)  ;
/* Return the length of the given multibyte character,
   putting its `wchar_t' representation in *PWC.  */
extern int mbtowc (wchar_t * __pwc,
		   const char * __s, size_t __n)  ;
/* Put the multibyte character represented
   by WCHAR in S, returning its length.  */
extern int wctomb (char *__s, wchar_t __wchar)  ;


/* Convert a multibyte string to a wide char string.  */
extern size_t mbstowcs (wchar_t *  __pwcs,
			const char * __s, size_t __n) ;
/* Convert a wide char string to multibyte string.  */
extern size_t wcstombs (char * __s,
			const wchar_t * __pwcs, size_t __n)
     ;



/* Determine whether the string value of RESPONSE matches the affirmation
   or negative response expression as specified by the LC_MESSAGES category
   in the program's current locale.  Returns 1 if affirmative, 0 if
   negative, and -1 if not matching.  */
extern int rpmatch (const char *__response)   ;






/* X/Open pseudo terminal handling.  */

/* Return a master pseudo-terminal handle.  */
extern int posix_openpt (int __oflag) ;



/* Put the 1 minute, 5 minute and 15 minute load averages into the first
   NELEM elements of LOADAVG.  Return the number written (never more than
   three, but may be less than NELEM), or -1 if an error occurred.  */
extern int getloadavg (double __loadavg[], int __nelem)
      ;


/* Define some macros helping to catch buffer overflows.  */




/* Copyright (C) 1991, 1996, 1997 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */

/*
 * SVID
 */


/* features.h Copyright FUJITSU LIMITED 2008-2010 */



/* Copyright (C) 1991-1993, 1995-2004, 2007 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */

/*
 *	ISO C99 Standard: 7.21 String handling	<string.h>
 */


/* features.h Copyright FUJITSU LIMITED 2008-2010 */




/* Get size_t and NULL from <stddef.h>.  */
/* stddef.h COPYRIGHT FUJITSU LIMITED 1998-2009                          */


 

 
 
 







/* Copy N bytes of SRC to DEST.  */
extern void *memcpy (void * __dest,
		     const void * __src, size_t __n)
      ;
/* Copy N bytes of SRC to DEST, guaranteeing
   correct behavior for overlapping strings.  */
extern void *memmove (void *__dest, const void *__src, size_t __n)
      ;


/* Copy no more than N bytes of SRC to DEST, stopping when C is found.
   Return the position in DEST one byte past where C was copied,
   or NULL if C was not found in the first N bytes of SRC.  */
extern void *memccpy (void * __dest, const void * __src,
		      int __c, size_t __n)
      ;



/* Set N bytes of S to C.  */
extern void *memset (void *__s, int __c, size_t __n)  ;

/* Compare N bytes of S1 and S2.  */
extern int memcmp (const void *__s1, const void *__s2, size_t __n)
       ;

/* Search N bytes of S for C.  */
extern void *memchr (const void *__s, int __c, size_t __n)
        ;





/* Copy SRC to DEST.  */
extern char *strcpy (char * __dest, const char * __src)
      ;
/* Copy no more than N characters of SRC to DEST.  */
extern char *strncpy (char * __dest,
		      const char * __src, size_t __n)
      ;

/* Append SRC onto DEST.  */
extern char *strcat (char * __dest, const char * __src)
      ;
/* Append no more than N characters from SRC onto DEST.  */
extern char *strncat (char * __dest, const char * __src,
		      size_t __n)  ;

/* Compare S1 and S2.  */
extern int strcmp (const char *__s1, const char *__s2)
       ;
/* Compare N characters of S1 and S2.  */
extern int strncmp (const char *__s1, const char *__s2, size_t __n)
       ;

/* Compare the collated forms of S1 and S2.  */
extern int strcoll (const char *__s1, const char *__s2)
       ;
/* Put a transformation of SRC into no more than N bytes of DEST.  */
extern size_t strxfrm (char * __dest,
		       const char * __src, size_t __n)
      ;



/* Duplicate S, returning an identical malloc'd string.  */
extern char *strdup (const char *__s)
       ;

/* Return a malloc'd copy of at most N bytes of STRING.  The
   resultant string is terminated even if no null terminator
   appears before STRING[N].  */



/* Find the first occurrence of C in S.  */
extern char *strchr (const char *__s, int __c)
       ;
/* Find the last occurrence of C in S.  */
extern char *strrchr (const char *__s, int __c)
       ;




/* Return the length of the initial segment of S which
   consists entirely of characters not in REJECT.  */
extern size_t strcspn (const char *__s, const char *__reject)
       ;
/* Return the length of the initial segment of S which
   consists entirely of characters in ACCEPT.  */
extern size_t strspn (const char *__s, const char *__accept)
       ;
/* Find the first occurrence in S of any character in ACCEPT.  */
extern char *strpbrk (const char *__s, const char *__accept)
       ;
/* Find the first occurrence of NEEDLE in HAYSTACK.  */
extern char *strstr (const char *__haystack, const char *__needle)
       ;


/* Divide S into tokens separated by characters in DELIM.  */
extern char *strtok (char * __s, const char * __delim)
      ;


/* Divide S into tokens separated by characters in DELIM.  Information
   passed between calls are stored in SAVE_PTR.  */
extern char *__strtok_r (char * __s,
			 const char * __delim,
			 char ** __save_ptr)
      ;
extern char *strtok_r (char * __s, const char * __delim,
		       char ** __save_ptr)
      ;





/* Return the length of S.  */
extern size_t strlen (const char *__s)
       ;





/* Return a string describing the meaning of the `errno' code in ERRNUM.  */
extern char *strerror (int __errnum) ;

/* Reentrant version of `strerror'.
   There are 2 flavors of `strerror_r', GNU which returns the string
   and may or may not use the supplied temporary buffer and POSIX one
   which fills the string into the buffer.
   To use the POSIX version, -D_XOPEN_SOURCE=600 or -D_POSIX_C_SOURCE=200112L
   without -D_GNU_SOURCE is needed, otherwise the GNU version is
   preferred.  */
/* Fill BUF with a string describing the meaning of the `errno' code in
   ERRNUM.  */
extern int __xpg_strerror_r (int __errnum, char *__buf, size_t __buflen)
      ;



/* We define this function always since `bzero' is sometimes needed when
   the namespace rules does not allow this.  */
extern void __bzero (void *__s, size_t __n)  ;

/* Copy N bytes of SRC to DEST (like memmove, but args reversed).  */
extern void bcopy (const void *__src, void *__dest, size_t __n)
      ;

/* Set N bytes of S to 0.  */
extern void bzero (void *__s, size_t __n)  ;

/* Compare N bytes of S1 and S2 (same as memcmp).  */
extern int bcmp (const void *__s1, const void *__s2, size_t __n)
       ;

/* Find the first occurrence of C in S (same as strchr).  */
extern char *index (const char *__s, int __c)
       ;

/* Find the last occurrence of C in S (same as strrchr).  */
extern char *rindex (const char *__s, int __c)
       ;

/* Return the position of the first bit set in I, or 0 if none are set.
   The least-significant bit is position 1, the most-significant 32.  */
extern int ffs (int __i)  ;

/* The following two functions are non-standard but necessary for non-32 bit
   platforms.  */

/* Compare S1 and S2, ignoring case.  */
extern int strcasecmp (const char *__s1, const char *__s2)
       ;

/* Compare no more than N chars of S1 and S2, ignoring case.  */
extern int strncasecmp (const char *__s1, const char *__s2, size_t __n)
       ;


/* Return the next DELIM-delimited token from *STRINGP,
   terminating it with a '\0', and update *STRINGP to point past it.  */
extern char *strsep (char ** __stringp,
		     const char * __delim)
      ;








/***************************************************************
* get bit template 
***************************************************************/
unsigned char
get_bit_tmpl(int ib){
	static unsigned char 
		bit_tmpl[]={
			0x80,0x40,0x20,0x10,
			0x08,0x04,0x02,0x01
		};
	if(ib>=0&&ib<8) return bit_tmpl[ib];
	else            return 0xff;
}

/***************************************************************
* get bit mask template
***************************************************************/
unsigned char
get_msk_tmpl(int ib){
	static unsigned char 
		msk_tmpl[]={
			0x7f,0xbf,0xdf,0xef,
			0xf7,0xfb,0xfd,0xfe
		};
	if(ib>=0&&ib<8) return msk_tmpl[ib];
	else            return 0xff;
}

/***************************************************************
* read bit data
***************************************************************/
int
read_bit(void *pc, int ib){
	unsigned char	tmpl,*pcpt;
	int				ibyt,ibit;

	ibyt=ib/8; ibit=ib%8;
	tmpl=get_bit_tmpl(ibit);
	pcpt=(unsigned char *)pc+ibyt;
	if(((*pcpt)&tmpl)==tmpl) return 1;
	else                     return 0;
}

/***************************************************************
* set 1 in bit data
***************************************************************/
void
set_bit1(void *pc, int ib){
	unsigned char	tmpl,*pcpt;
	int				ibyt,ibit;

	ibyt=ib/8; ibit=ib%8;
	tmpl=get_bit_tmpl(ibit);
	pcpt=(unsigned char *)pc+ibyt;
	(*pcpt)|=tmpl;
}

/***************************************************************
* set 0 in bit data
***************************************************************/
void
set_bit0(void *pc, int ib){
	unsigned char	tmpl,*pcpt;
	int				ibyt,ibit;

	ibyt=ib/8; ibit=ib%8;
	tmpl=get_msk_tmpl(ibit);
	pcpt=(unsigned char *)pc+ibyt;
	(*pcpt)&=tmpl;
}

/***************************************************************
* MOVEB for C
***************************************************************/
void
moveb_c(void* dst, int idst, void *src, int isrc, int leng){
	int	ic;
	for(ic=0;ic<leng;ic++){
		if(read_bit(src,isrc+ic)==1) set_bit1(dst,idst+ic);
		else                         set_bit0(dst,idst+ic);
	}
}

/***************************************************************
* MOVEB for Fortran
***************************************************************/
void
MOVEB
(void* dst, int *idst, void *src, int *isrc, int *leng){
	moveb_c(dst,*idst-1,src,*isrc-1,*leng);
}

/***************************************************************
* MOVEC for C
***************************************************************/
void
movec_c(void* dst, int idst, void *src, int isrc, int leng){
	if(leng>0){
		unsigned char	*cdst,*csrc;
		cdst=(unsigned char *)dst+idst;
		csrc=(unsigned char *)src+isrc;
		memcpy(cdst,csrc,leng);
	}
}

/***************************************************************
* MOVEC for Fortran
***************************************************************/
void
MOVEC
movec_
(void* dst, int *idst, void *src, int *isrc, int *leng)
{
	movec_c(dst, *idst-1, src, *isrc-1, *leng);
}

/***************************************************************
* show bit data for test
***************************************************************/
void show_bit(void *pc,int leng){
	int ic;
	for(ic=0;ic<leng;ic++){
		printf("%i",read_bit(pc,ic));
		if(ic%4==3) printf(" ");
	}
	printf("\n");
}
/***************************************************************
* main routien foe test
***************************************************************/
/*
void
main(void){
	int	src,dst,ic;
	src=65535;
	dst=0;
	show_bit(&src,16);
	show_bit(&dst,16);
	moveb((void *)&dst,4,(void *)&src,4,8);
	show_bit(&dst,16);
}
*/
