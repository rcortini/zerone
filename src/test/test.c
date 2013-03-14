#include <glib.h>
#include <stdio.h>
#include <string.h>

typedef struct {
   trienode *trie;
   hitlist  *hits;
} fixture;

void
setup(
   fixture *f,
   gconstpointer test_data
)
{
   return;
}

void
test_search(
   fixture *f,
   gconstpointer ignore
)
{
   search(f->trie, "AAA", 0, f->hits);
   g_assert(f->hits->n_hits == 1);

   clear_hitlist(f->hits);
   search(f->trie, "AAA", 1, f->hits);
   g_assert(f->hits->n_hits == 3);

   clear_hitlist(f->hits);
   search(f->trie, "AAA", 2, f->hits);
   g_assert(f->hits->n_hits == 4);
}


void
test_seq(
   fixture *f,
   gconstpointer ignore
)
{
   char buffer[8];
   trienode *node = insert(f->trie, "AAATGC", 0);
   g_assert(strcmp(seq(node, buffer, 8), "AAATGC") == 0);

   node = insert(f->trie, "GYATC", 0);
   g_assert(strcmp(seq(node, buffer, 8), "GNATC") == 0);

   search(f->trie, "AAA", 0, f->hits);
   g_assert(strcmp(seq(f->hits->node[0], buffer, 8), "AAA") == 0);
}


void
teardown(
   fixture *f,
   gconstpointer test_data
)
{
   destroy_trie(f->trie);
   destroy_hitlist(f->hits);
}


void
test_run
(void)
{
   FILE *inputf = fopen("input_test_file", "r");
   FILE *outputf = fopen("/dev/null", "w");
   starcode(inputf, outputf, 0);
}



int
main(
   int argc,
   char **argv
)
{
   g_test_init(&argc, &argv, NULL);
   g_test_add("/trie/search", fixture, NULL, setup, test_search, teardown);
   g_test_add("/trie/seq", fixture, NULL, setup, test_seq, teardown);
   g_test_add_func("/starcode/run", test_run);
   return g_test_run();
}
