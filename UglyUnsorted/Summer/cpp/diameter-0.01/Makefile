FILES = alloc.c \
        pick.c \
        print.c \
        rand.c \
        read.c \
        util.c \
        exact-diam.c \
        apprx-diam.c \
	peled-like.c 


all : gener-points diameter-test

CC=cc
OPT=-O

BINDIR = .

clean :
	rm -f ${BINDIR}/gener-points ${BINDIR}/diameter-test
 
gener-points : ${FILES} gener-points.c
	${CC} ${OPT}  -I. -o ${BINDIR}/gener-points gener-points.c ${FILES} -lm
	rm -f *.o

diameter-test : ${FILES} diameter-test.c 
	${CC} ${OPT}  -I. -o ${BINDIR}/diameter-test diameter-test.c ${FILES} -lm
	rm -f *.o

