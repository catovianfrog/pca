char *p;
int c;
enum states { DULL, IN_WORD, IN_STRING } state = DULL;

for (p = buffer; *p != '\0'; p++) {
    c = (unsigned char) *p; /* convert to unsigned char for is* functions */
    switch (state) {
    case DULL: /* not in a word, not in a double quoted string */
        if (isspace(c)) {
            /* still not in a word, so ignore this char */
            continue;
        }
        /* not a space -- if it's a double quote we go to IN_STRING, else to IN_WORD */
        if (c == '"') {
            state = IN_STRING;
            start_of_word = p + 1; /* word starts at *next* char, not this one */
            continue;
        }
        state = IN_WORD;
        start_of_word = p; /* word starts here */
        continue;

    case IN_STRING:
        /* we're in a double quoted string, so keep going until we hit a close " */
        if (c == '"') {
            /* word goes from start_of_word to p-1 */
            ... do something with the word ...
            state = DULL; /* back to "not in word, not in string" state */
        }
        continue; /* either still IN_STRING or we handled the end above */

    case IN_WORD:
        /* we're in a word, so keep going until we get to a space */
        if (isspace(c)) {
            /* word goes from start_of_word to p-1 */
            ... do something with the word ...
            state = DULL; /* back to "not in word, not in string" state */
        }
        continue; /* either still IN_WORD or we handled the end above */
    }
}
