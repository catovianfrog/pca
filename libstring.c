/**********************************************************************
* LIBSTRING
*
* String manipulation library
* v0.1	June 2014
* (C) Bruno Charrière, see http://www.gnu.org/licenses/gpl.html
**********************************************************************/
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include "libstring.h"


//---------- Private routine used by tokenize() --------------------------------
char** get_token(int word_beg, int word_end, int *word_count, char *buffer, char **words) {
    (*word_count)++;
    buffer[word_end+1]='\0'; // terminate string
    words=realloc(words,*(word_count)*sizeof(words)); //extend size of word array
    words[*word_count-1]=calloc(1,word_end-word_beg+1);
    strcpy(words[*word_count-1],&buffer[word_beg]);	// token points to strt of word
    return words;
}
 /**********************************************************************
*	    Tokenize
* This function splits a string into words or substrings
**********************************************************************/
char** tokenize(char *buf,int *word_count)
{
    enum    states {between_words, in_word, in_quote, in_dquote} state;
    int	    ch,pos,word_beg;
    char    **tokens;
    char*   err_msg;

    tokens=NULL;
    state=between_words;
    *word_count=0;
    for(pos=0; buf[pos]!='\0'; pos++) {	//parse the buffer string
	ch=buf[pos];
	    switch (state) {
    	    case between_words:
		if(isspace(ch)) continue;// still not in a word, ignore this char
		// ch is not a space. Enter either a quoted string or a word
                if(ch == '\'') {
		    state=in_quote;
		    word_beg=pos+1; //word starts after the quote
		    continue;
		}
                if(ch == '"') {
		    state=in_dquote;
		    word_beg=pos+1; //word starts after the quote
		    continue;
		} // not a space nor a quote: in a word
		state=in_word;
		word_beg=pos;	// word starts here
		continue;
	    case in_quote:	// keep going until next quote
		if(ch == '\'') {
		    tokens=get_token(word_beg, pos-1, word_count, buf, tokens);
		    state=between_words;
		}
		continue; //either still in string or between words
	    case in_dquote:	// keep going until next quote
		if(ch == '"') {
		    tokens=get_token(word_beg, pos-1, word_count, buf, tokens);
		    state=between_words;
		}
		continue;	//either still in string or between words
	    case in_word:	//keep going until next space
		if(isspace(ch)) {
		    tokens=get_token(word_beg, pos-1, word_count, buf, tokens);
		    state=between_words;
		}
		if(ch=='\'' || ch=='\"') {
		    err_msg="String badly constructed. Quote within a word.";
		    goto error;
		}
		continue;
	    }
    }
    // end of buffer is reached ('\0');
    if(state == in_word) tokens=get_token(word_beg, pos-1, word_count, buf, tokens);
    if(state == in_quote || state== in_dquote) {
        err_msg="String badly constructed. Uneven number of quotes.";
	goto error;
    }
    return tokens;

error:
    fprintf(stderr, "******** ERROR: %s\n",err_msg);
    *word_count=0;
    tokens=NULL;
    return tokens;
}

