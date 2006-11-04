#include <stdio.h>
#include <readline/readline.h>
#include <readline/history.h>

char *pbc_getline(void)
{
    char *line = readline(NULL);
    if (line && *line) add_history(line);
    return line;
}
