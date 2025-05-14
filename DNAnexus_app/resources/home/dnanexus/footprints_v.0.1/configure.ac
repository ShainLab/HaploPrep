# initial information about the project
AC_INIT([foots], [0.1], [eleanorrose.fewings@ucsf.edu])
AM_INIT_AUTOMAKE

# Check for programmes
AC_PROG_CC
AC_PATH_PROG(samv019, samtools --version, no)
AC_SUBST(samv019, "$(dirname $samv019)")
if test "$samv019" = "no"; then
   AC_MSG_ERROR([samtools-0.1.19 not found on PATH!])
else
   AC_MSG_RESULT([samtools-0.1.19 path: $samv019])
fi

# Set configure files
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
