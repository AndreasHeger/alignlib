dnl #########################################################################
dnl Check command line option --enable-DEBUG
AC_DEFUN([adl_ENABLE_DEBUG],
[
AC_PROVIDE([$0])
 
ac_cv_enable_debug=no
 
dnl command line options and --help descriptions
AC_ARG_ENABLE( debug,
[  --enable-debug		compile with option DEBUG \[default=no\] ],
ac_cv_enable_debug=$enableval,	
ac_cv_enable_debug=no	
)

if test x"$ac_cv_enable_debug" = xyes; then
	CXXFLAGS="$CXXFLAGS -DDEBUG"
	AC_SUBST(CXXFLAGS)

	AC_MSG_CHECKING("Enable option DEBUG")
	AC_MSG_RESULT(yes)
fi
]
)

dnl Handle the --enable-html-doc option.
dnl
dnl The following snippet goes into doc/Makefile.am
dnl ------------------------
dnl CLEANFILES = yourfile.html
dnl
dnl if MAKE_HTML
dnl hdir = @htmldir@
dnl h_DATA = yourfile.html
dnl endif
dnl
dnl SUFFIXES = .html
dnl
dnl .texi.html:
dnl         @cd $(srcdir) && rm -f $@ $@-[0-9] $@-[0-9][0-9]
dnl         cd $(srcdir) && $(MAKEINFO) --html `echo $< | sed 's,.*/,,'`
dnl ------------------------
dnl
dnl Written by Alexandre Duret-Lutz <duret_g@epita.fr>
 
AC_DEFUN([adl_ENABLE_HTML_DOC],
[
AC_ARG_ENABLE([html-doc],
[  --enable-html-doc=DIR       build and install html documentation (install via make install-data)])
  if test "${enable_html_doc-no}" != no; then
    if test "x$enable_html_doc" = xyes; then
      htmldir="$prefix/doc/$PACKAGE"
    else
      htmldir="$enable_html_doc"
    fi
    AC_SUBST([htmldir])
  else
    htmldir='\<none\>'
  fi
  AM_CONDITIONAL(MAKE_HTML, [test "${enable_html_doc-no}" != no])
  ]
  )     







