CFLAGS=-Wall -g -DNDBUG -O0
LDFLAGS = -lm

pca: 	pca.c libal.o dbg.h libstring.o
	gcc $(CFLAGS)  -c pca.c  $(LDFLAGS)   
	gcc $(CFLAGS)  pca.o libal.o libstring.o -o pca  $(LDFLAGS)   
	markdown README.md >readme.html

libal.o: libal.c dbg.h
	gcc $(CFLAGS)  -c libal.c  $(LDFLAGS)   

libstring.o:	libstring.c
	gcc $(CFLAGS)  -c libstring.c    

clean:
	rm -f *.o
	rm -f test
	rm -f pca
	rm -f al
	rm -f tmp.*
	rm -f temp.*   
	rm -f *~

	

