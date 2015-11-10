EXEC_NAME=fox
CC=mpiCC

# CFLAGS= -Wall -Wno-write-strings -O0 -g
CFLAGS= -w -Wno-write-strings -O3 -Wall -g

SRC =                   \
	main.cpp	\
	matrix.cpp	\
	fox.cpp		\
	floyd.cpp	\
	seqfox.cpp	\
	seqfloyd.cpp

OBJ =  ${SRC:.cpp=.o}

#------------------------------------------------------------

all: ${EXEC_NAME}

${EXEC_NAME}: ${OBJ}
	${CC} ${CFLAGS} ${CLIBS} -o ${EXEC_NAME} ${OBJ}

%.o: %.cpp
	${CC} ${CFLAGS} -c -o $@ $+

clean:
	rm ${EXEC_NAME} *.o *~ *# -rf
