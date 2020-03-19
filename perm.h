#ifndef PERM_H

    #define PERM_H

    typedef struct perm_node_struct *Node;
    
    typedef struct perm_node_struct {
    
        int tag;
        Node next;
        Node last;
    
    } perm_node;
    
    typedef struct perm_list_struct {
        
        int **L; /* list */
        int n; /* number of rows in list */
        int k; /* number of elements in each bin */
        int m; /* number of bins */
        
    } perm_list;

    void perm_get(perm_list *, perm_node *);
    void perm_put(perm_list *, perm_node *, perm_node *, int, int);
    void perm_pop(perm_list *, perm_node *);
    int perm_power(int, int);
    int **perm(int, int);
    
#endif

