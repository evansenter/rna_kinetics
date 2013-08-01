/* From libRNA.a */
float energy_of_struct(char *, char *);

void energy_of_struct_r(char **sequence, char **structure, double *energy) {
  *energy = energy_of_struct(*sequence, *structure);
}

