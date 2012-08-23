#ifndef EFPMD_MD_H
#define EFPMD_MD_H

struct efp;
struct config;
struct md;

struct md *md_create(struct efp *, const struct config *);
void md_update_step_nve(struct md *);
void md_update_step_nvt(struct md *);
void md_print_info(struct md *);
void md_print_geometry(struct md *);
void md_shutdown(struct md *);

#endif /* EFPMD_MD_H */
