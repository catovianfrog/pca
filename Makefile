CFLAGS=-Wall -g -DNDBUG -O0
LDFLAGS = -lm

pca: 	pca.c libal.c
	gcc $(CFLAGS)  pca.c -o pca $(LDFLAGS)   

libal: 	libal.c dbg.h
	gcc $(CFLAGS)  -c libal.c  $(LDFLAGS)   


clean:
	rm -f *.o
	rm -f test
	rm -f pca
	rm -f al
	

