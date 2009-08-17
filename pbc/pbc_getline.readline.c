#include <stdio.h>
#include <readline/readline.h>
#include <readline/history.h>

char *pbc_getline(const char *prompt)
{
    char *line = readline(prompt);
    if (line && *line) add_history(line);
    return line;
}
