#ifndef EMFIELD_H
#define EMFIELD_H


/** structure with field information */
struct EMfield {
    // field arrays: 4D arrays
    
    /* Electric field defined on nodes: last index is component */
    FPfield*** Ex;
    FPfield*** Ey;
    FPfield*** Ez;
    /* Magnetic field defined on nodes: last index is component */
    FPfield*** Bxn;
    FPfield*** Byn;
    FPfield*** Bzn;
    
    
};

/** allocate electric and magnetic field */
inline void field_allocate(struct grid* grd, struct EMfield* field){
    
    // E on nodes
    field->Ex  = newArr3(FPfield, grd->nxn, grd->nyn, grd->nzn);
    field->Ey  = newArr3(FPfield, grd->nxn, grd->nyn, grd->nzn);
    field->Ez = newArr3(FPfield, grd->nxn, grd->nyn, grd->nzn);
    // B on nodes
    field->Bxn = newArr3(FPfield, grd->nxn, grd->nyn, grd->nzn);
    field->Byn = newArr3(FPfield, grd->nxn, grd->nyn, grd->nzn);
    field->Bzn = newArr3(FPfield, grd->nxn, grd->nyn, grd->nzn);
    
}

/** deallocate electric and magnetic field */
inline void field_deallocate(struct grid* grd, struct EMfield* field){
    // E deallocate 3D arrays
    delArr3(field->Ex, grd->nxn, grd->nyn);
    delArr3(field->Ey, grd->nxn, grd->nyn);
    delArr3(field->Ez, grd->nxn, grd->nyn);
    // B deallocate 3D arrays
    delArr3(field->Bxn, grd->nxn, grd->nyn);
    delArr3(field->Byn, grd->nxn, grd->nyn);
    delArr3(field->Bzn, grd->nxn, grd->nyn);
}

#endif
