#ifndef WORKQUEUE_H
#define WORKQUEUE_H

#include "jurassic.h"

/* Module containing all variables and functions related to a work-queue for scattering */

#define Queue_Inactive -1
#define Queue_Prepare   2 
#define Queue_Execute   8
#define Queue_Collect  32

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

int init_queue(queue_t *q, int size);
int push_queue(queue_t *q, void* inp, void* out, int ir);
int get_queue_item(queue_t *q, void** inp, void **out, int *ir, int index);
int pop_queue(queue_t *q, void** inp, void **out, int *ir);

#endif
