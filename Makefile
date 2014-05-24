CFLAGS=-Wall -g -DNDBUG -O0
LDFLAGS = -lm

pca: 	pca.c
	gcc $(CFLAGS)  pca.c -o pca $(LDFLAGS)   

al: 	al.c
	gcc $(CFLAGS)  al.c -o al $(LDFLAGS)   

test:   test.c
	cc test.c -o test -lm

clean:
	rm -f *.o
	rm -f qr
	rm -f mat
	rm -f test
	rm -f pca
	

