CFLAGS=-Wall -g -DNDBUG -O0
LDFLAGS = -lm

pca: 	pca.c libal.o dbg.h
	gcc $(CFLAGS)  -c pca.c  $(LDFLAGS)   
	gcc $(CFLAGS)  pca.o libal.o -o pca  $(LDFLAGS)   

libal.o: libal.c dbg.h
	gcc $(CFLAGS)  -c libal.c  $(LDFLAGS)   


clean:
	rm -f *.o
	rm -f test
	rm -f pca
	rm -f al
	rm -f tmp.*
	rm -f temp.*   
	rm -f *~

	

