SUBDIR = common src

.if !defined(WITHOUT_EFPMD)
SUBDIR += efpmd
.endif

.if !defined(WITHOUT_TESTS)
SUBDIR += tests
.endif

.include <subdir.mk>
