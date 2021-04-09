#define dgsand_calloc(t,n) (t *) calloc((n),sizeof(t));
#define dgsand_alloc(t,n) (t *)malloc((n)*sizeof(t))
#define dgsand_free(t) free(t)
