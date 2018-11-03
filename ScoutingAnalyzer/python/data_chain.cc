TChain *data_chain(const char *flist, const char *bname) {
    FILE *z = fopen(flist, "r");
   
    TChain *dc = new TChain(bname);    
    char *line = nullptr;
    size_t line_c = 0;
    ssize_t len = 0;
    
    while((len = getline(&line, &line_c, z)) > 0) {
                line[len - 1] = '\0';
                dc->Add(line);
        }

   return dc;
}
