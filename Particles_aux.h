#ifndef PARTICLES_AUX_H
#define PARTICLES_AUX_H

struct particles_aux {
    
    /** species ID: 0, 1, 2 , ... */
    int species_ID;
    
    /** maximum number of particles of this species on this domain. used for memory allocation */
    long npmax;
    /** number of particles of this species on this domain */
    long nop;
    
    /** densities carried nop,2,2,2*/
    FPpart (*rho_p)[2][2][2];
    FPpart (*Jx)[2][2][2]; FPpart (*Jy)[2][2][2]; FPpart (*Jz)[2][2][2];
    FPpart (*pxx)[2][2][2]; FPpart (*pxy)[2][2][2]; FPpart (*pxz)[2][2][2];
    FPpart (*pyy)[2][2][2]; FPpart (*pyz)[2][2][2]; FPpart (*pzz)[2][2][2];
    
    
    /** cell index: ix, iy, iz */
    int* ix_p; int* iy_p; int* iz_p;
    
};

/** allocate particle arrays */
inline void particle_aux_allocate(struct particles* part, struct particles_aux* part_aux, int is){
    
    // set species ID
    part_aux->species_ID = is;
    // number of particles
    part_aux->nop = part->nop;
    // maximum number of particles
    part_aux->npmax = part->npmax;
    
    long npmax = part->npmax;
    
    // allocate densities brought by each particle
    part_aux->rho_p  = new FPpart[part->npmax][2][2][2];
    part_aux->Jx  = new FPpart[part->npmax][2][2][2];
    part_aux->Jy  = new FPpart[part->npmax][2][2][2];
    part_aux->Jz  = new FPpart[part->npmax][2][2][2];
    part_aux->pxx  = new FPpart[part->npmax][2][2][2];
    part_aux->pxy  = new FPpart[part->npmax][2][2][2];
    part_aux->pxz  = new FPpart[part->npmax][2][2][2];
    part_aux->pyy  = new FPpart[part->npmax][2][2][2];
    part_aux->pyz  = new FPpart[part->npmax][2][2][2];
    part_aux->pzz  = new FPpart[part->npmax][2][2][2];
    
    // cell index
    part_aux->ix_p = new int[part->npmax];
    part_aux->iy_p = new int[part->npmax];
    part_aux->iz_p = new int[part->npmax];
    
    
}

inline void particle_aux_deallocate(struct particles_aux* part_aux){
    // deallocate auxiliary particle variables needed for particle interpolation
    delete [] part_aux->rho_p;
    delete [] part_aux->Jx;
    delete [] part_aux->Jy;
    delete [] part_aux->Jz;
    delete [] part_aux->pxx;
    delete [] part_aux->pxy;
    delete [] part_aux->pxz;
    delete [] part_aux->pyy;
    delete [] part_aux->pyz;
    delete [] part_aux->pzz;
    
    
}

#endif

