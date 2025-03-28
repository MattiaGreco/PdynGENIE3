/* Wrapper of Pierre Geurts' C code for learning trees
 *
 * To compile:
 *   R CMD SHLIB PdynGENIE3.c
 *   or
 *   R CMD SHLIB -g PdynGENIE3.c      (This will allow you to debug the shared 
 *                                     library via for instance gdb or valgrind.)
 *
 *
 * Modifications by Van Anh:
 *
 * August 2016:
 *      - Merged everything into one file
 *      - Removed functions that are not necessary for the R package
 *      - Made the following changes (otherwise it generates a warning when doing the R CMD check):
 *             - removed all the fprintf
 *             - modified a part of the code so that there is no exit(0) anymore (a flag_continue is used instead)
 *             - replaced rand() with unif_rand()
 *
 * Modifications by Mattia Greco and Olivier Martin:
 *
 * July 2024:
 *      - A given tree now considers only the same K features for all splits
 *        Cf. find_the_best_split_among_k_preselected  
 *        One can return to the previous method by using find_the_best_split_among_k instead
 *      - We choose the K features according to priors that are passed
 *      - Additional properties of the forest are computed like its prediction
 */

#include <R.h>
#include <Rinternals.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/*#include <string.h>*/
/*#include <ctype.h>*/


/**********************************************/
/* From here: what used to be in tree-model.h */
/**********************************************/


/*
 * Author: Pierre Geurts (p.geurts@ulg.ac.be)
 * (c) 2002-2010
 *
 * (Please do not distribute without permission of the author)
 *
 */

/* MACRO, CONSTANTES */

/* precision des differentes fonction de calcul */

#define MAX_OPEN_NODES 15000

#define SCORE_TYPE double
/* permet de definir la core_table en char pour les grosses BDs */
#define CORETABLE_TYPE float
/* #define CORETABLE_TYPE unsigned char */

#define MAX_NUMBER_OF_ATTRIBUTES 300000

/* pour les variables symboliques */

#define MAX_NUMBER_OF_SYMBOLIC_VALUES 1024      /* 32 * 32 */
#define MAX_NUMBER_OF_SYMBOLIC_VALUES_DIV_32 32 /* C'est à dire maximum 128 valeurs */

#define BITN(x,n) (x.i[n/32]>>(n%32))%2
#define SET_BITN(x,n) x.i[n/32]|=(1<<(n%32))

#define SYMBOLICp(att) attribute_descriptors[att]>0
#define NUMERICALp(att) attribute_descriptors[att]==0
#define NB_VAL(att) attribute_descriptors[att]

#define MAX(x,y) (x<y? y: x)
#define MIN(x,y) (x<y? x: y)

/* redefinition du malloc (pour l'interface matlab) */
#define MyMalloc malloc
#define MyFree free

/* definition de types */

union threshold_type {
  unsigned int i[MAX_NUMBER_OF_SYMBOLIC_VALUES_DIV_32];
  float f;
};

/* prototype de fonctions */

/* fonction qui calcule des stats pour le vecteur et qui place le resultat dans table_score */
/* Specifically, initializes the total weight of the data points (unchanged when bagging because data are weighted)  */
/* as well as Sum_i weight_i*y_i and Sum_i weight_i*y_i**2 in the first row of the table_score */
void (*summarize_vector)(int *vector, int start, int end);
/* function that returns an index to get prediction for a leaf and adds to prediction_values array */
int (*make_leaf_prediction)();  /* usage is   prediction_values[prediction[node_in_tree_in_forest]][0] for 0th prediction */
                                /* the predicted value is the average of the y data at that leaf */
/* fonction qui renvoie true si on doit arreter la construction */
int (*stop_splitting_criterio)();  /* Only called for deciding to split the root node */
/* fonction qui renvoie true si le test choisi n'est pas significatif */
int (*not_significant_test)();
/* fonction qui calcule le score a partir d'une table de score */
SCORE_TYPE (*compute_score_from_table)();
/* une fonction qui recherche un test pour splitter le noeud */
void (*find_a_split)(int *ls_vector, int start, int end);
/* une fonction qui recherche un seuil pour un attribut numerique */
void (*find_a_threshold_num)(int att, int *ls_vector, int start, int end);
/* une fonction qui recherche un seuil pour un attribut symbolique */
void (*find_a_threshold_symb)(int att, int *ls_vector, int start, int end);

float *(*get_tree_prediction_vector)(int tree, int obj);  /* fetch prediction(s) for a given tree and obj */
                                                          /* assuming prediction_values array was saved   */


void write_one_tree(int tree, FILE *fp);
/* void ocm_write_one_tree(int tree, FILE *fp); */

/* Sorting function */

void quicksort_ls_vector(int ls_vector[], int start, int end, int att);
int separate_ls_vector_local(int best_attribute, union threshold_type best_threshold, int *ls_vector, int start, int end);

/* generique */
int build_one_tree();
void find_a_split_at_random(int *ls_vector, int start, int end);
void find_a_threshold(int att, int *ls_vector, int start, int end);
int check_test(int att, float val, union threshold_type threshold);

/* utilitaires */
int get_random_integer(int max_val);
float get_random_float();

float getattval(int obj, int att);

/* interface */
void clean_all_trees();
void clean_all_ensemble();

float make_ls_vector_bagging(int tree);
float make_ls_vector_identity(int tree);

int get_tree_nb_nodes(int tree);

void free_multiregr_table_score();





/***********************************************/
/* From here: what used to be in tree-model.c  */
/***********************************************/



/*
 * Author: Pierre Geurts (p.geurts@ulg.ac.be)
 * (c) 2002-2010
 *
 * (Please do not distribute without permission of the author)
 *
 * Cleaning by Vân Anh:
 *
 * July 2016:
 *		- Removed best first option
 *      - Removed significance test (f-test) at each node
 */


void init_threshold_type(union threshold_type *u) {
  int i;
  for (i=0;i<MAX_NUMBER_OF_SYMBOLIC_VALUES_DIV_32; i++)
    (*u).i[i]=0;
}

void add1_threshold_type(union threshold_type *u) {
  int i=0;
  do {
    (*u).i[i]++;
  } while (((*u).i[i++]==0)&&(i<=MAX_NUMBER_OF_SYMBOLIC_VALUES_DIV_32));
}

/**********************/
/* VARIABLES GLOBALES */
/**********************/

/* DEFINITION DE LA "TABLE" POUR STOCKER LES ARBRES */

int size_current_tree_table=0;
int size_current_tree_table_pred=0;

union threshold_type *threshold=NULL; /* le seuil associe a un noeud (sur 64 BITS pour reprÈsenter plus de valeurs
				                      * symboliques */
int *left_successor=NULL; /* la position RELATIVE du successeur par rapport */
int *right_successor=NULL; /* a ce noeud */
int *tested_attribute=NULL; /* l'attribut teste en ce noeud */
int *prediction=NULL; /* l'index de la prediction associe au noeud */
float *node_weight=NULL; /* un poids pour le test en dual perturb and combine */
/* Pour missing values: la taille de la feuille en nb d'elements */
float *node_size=NULL;

float **prediction_values=NULL;
int index_prediction_values=-1;

int index_nodes=-1;

int stack_open_nodes[MAX_OPEN_NODES][3];

int index_stack_open_nodes=-1;

/* table des valeurs */

CORETABLE_TYPE *core_table=NULL;

int nb_obj_in_core_table;  /* total number of objects (samples) */
                           /* The first nb_obj_in_core_table - nb_obj_in_validation_table are to be used for constructing trees */
int nb_obj_in_validation_table;    /* number of objects (samples) to be used for validating trees, taken from bottom of data set */
/* vector with attribute descriptor (0 if numerical, the number of values
 * if symbolic)
 */
int *attribute_descriptors=NULL;
int length_attribute_descriptors;

/* parametres generaux pour la construction */

int goal; /* goal_classification ou goal_regression */

int *attribute_vector=NULL; /* indice des attributs */
int nb_attributes;


int rf_k=-1;  // Number of attributes considered in splits

float *priors=NULL;
float *cumulative_priors=NULL;
float *PRESELECTED_ATTRIBUTES=NULL;   // Contains the K attributes to use in all splits of a given tree
                                      // These are chosen via priors but need not be distinct

int *current_learning_set=NULL; /* l'ensemble d'apprentissage */
int current_learning_set_size;  /* Number of samples (input,output) vectors used after bagging */

int *object_mapping=NULL;
SCORE_TYPE *object_weight=NULL;
int global_learning_set_size=0; /* taille des table object_mapping et object_weight: peut Ítre plus grand
				   que current_learning_set_size si certains objets sont de poids nuls*/
                                /* Number of samples (input,output) vectors in total data  */

int min_node_size; /* nombre min d'objets dans un noeud */
int MAX_NUM_NODES_IN_TREE = -1;  /* If -1, no limit; 7 if 2 drivers that are monotonic, 11 if one driver is and the other is bell shaped, etc */

/* parametres du tri: *
 * par defaut, on trie chaque fois localement */
void (*sort_ls_vector)(int *ls_vector, int start, int end, int att)=quicksort_ls_vector;
int (*separate_ls_vector)(int best_attribute, union threshold_type best_threshold, int *ls_vector, int start, int end)=separate_ls_vector_local; 

/* methodes d'ensemble */

float (*make_ls_vector)(int tree);
int number_of_ensemble_terms=50;
int save_ensemble_while_growing=0;  /* 1 to write tree data to file, 0 not to */
int store_ensemble=1;


/* parametres generaux pour la classification */

int nb_classes; /* nombre de classes pour la classification */  /* Not used except for write_one_tree that has been deactivated in PdynGENIE3! */

/* parametres generaux pour la regression */ 

/* ... */

/* pour le calcul du score */

SCORE_TYPE **table_score=NULL; /* [3][MAX_NUMBER_OF_PREDICTION_VALUES+1]; */
SCORE_TYPE **table_score_symb=NULL; /* [MAX_NUMBER_OF_SYMBOLIC_VALUES][MAX_NUMBER_OF_PREDICTION_VALUES+1]; */
int nb_of_predictions;

/* on pourrait utiliser table_score et table_score_symb mais je prefere utiliser
 * des tables dediees pour que MAX_GOAL_MULTIREGR soit independant de MAX_NUMBER_OF_NODES
 */
SCORE_TYPE **table_score_multiregr=NULL; /* [3][2*MAX_GOAL_MULTIREGR+1]; */
SCORE_TYPE **table_score_symb_multiregr=NULL; /* [MAX_NUMBER_OF_SYMBOLIC_VALUES][2*MAX_GOAL_MULTIREGR+1]; */

/* parametre de stop splitting pour les arbres */

float h_min;

/* variables globales pour le calcul des meilleures score */

int best_attribute;
union threshold_type best_threshold;
SCORE_TYPE best_threshold_score;

union threshold_type current_threshold;
SCORE_TYPE current_threshold_score;

/* variables globales pour les stop splitting criterio a posteriori (not-significant) */
SCORE_TYPE current_threshold_info;
SCORE_TYPE best_threshold_info;
SCORE_TYPE total_weighted_variance_at_tree_root;
int keep_split;
float MIN_FRAC_OF_TOTAL_WEIGHTED_VARIANCE_IN_SPLIT=0.0;  // minimum improvement when doing a split (compared to total variance)
float MIN_FRAC_OF_VARIANCE_BEFORE_SPLIT=0.0;    // minimum improvement when doing a split (compared to previous variance)
SCORE_TYPE WEIGHTED_VARIANCE_EXPLAINED;                   // for each tree, gives the (weighted) variance explained on its (weighted) LS set

/* pour le calcul du score */
SCORE_TYPE info;
SCORE_TYPE v_tot;   /* computed in stop_splitting_criterio */
SCORE_TYPE v_min=0.0;

/* Global variables for effective number of attributes used each tree */

float *Effective_Num_Attributes_In_Trees;
float Neff;

/* Global variables for the square prediction error on validation data and the learning set size for each tree  */
int CURRENT_itree;   /* The number of the tree (0, 1, 2...) currently being treated */
float *SqE_on_validation; /* one entry for each tree, normalized by the variance of that validation data */
float *LS_size;

/* Global float for the square error of the tree-averaged predictions, averaged over validation iobj */
/* divided by the variance of that validation data */
float *pForest_Prediction_SS_on_validation_to_Return=NULL;  /* points to a single float */

float FOREST_QUALITY;



/***************/
/* CODE ARBRES */
/***************/

/* BUILD_ONE_TREE */
/*----------------*/
/* fonction generale de construction d'un arbre */

CORETABLE_TYPE *core_table_y;  // also defined elsewhere

int build_one_tree() {
  int tree;
  int nb_tests=0;
  int root_node = 0;
  float y_value, sum_y, sum_y2, weight, wtot, variance_full_data;

// DEBUG:
//  printf("\n In build_one_tree \n");
//  for(int ii = 0; ii < global_learning_set_size; ++ii) printf("y[%d] = %e  ",ii,core_table_y[ii]);
  sum_y=0.0;
  sum_y2=0.0;
  for(int ii = 0; ii < global_learning_set_size; ++ii) {
    y_value = core_table_y[ii];
    sum_y += y_value;
    sum_y2 += y_value*y_value;
  }
  variance_full_data = (sum_y2 - (sum_y*sum_y)/global_learning_set_size)/global_learning_set_size;
//  printf("###  total variance of full data set = %e  ",(sum_y2 - (sum_y*sum_y)/global_learning_set_size)/global_learning_set_size);

  sum_y=0.0;
  sum_y2=0.0;
  wtot=0.0;
  for(int ii = 0; ii < current_learning_set_size; ++ii) {
    weight = object_weight[current_learning_set[ii]];
    wtot += weight;
    y_value = core_table_y[current_learning_set[ii]];
    sum_y += weight*y_value;
    sum_y2 += weight*y_value*y_value;
  }
  total_weighted_variance_at_tree_root=(sum_y2 - (sum_y*sum_y)/wtot);
// printf("LS-based first row of table_score  = %e  %e  %e\n",wtot,sum_y,sum_y2);
// printf("total_weighted_variance_of_LS = %e  ####\n",total_weighted_variance_at_tree_root);

  WEIGHTED_VARIANCE_EXPLAINED = 0.0;

// If the tree is built using a fixed set of K attributes, construct the PRESELECTED_ATTRIBUTES array

  int random_att_pos;
  int lleft;
  int rright;
  float x_rand;

  for(int ith_preselected = 0; ith_preselected < rf_k; ++ith_preselected) {
    lleft = 0;
    rright = nb_attributes-1;
    x_rand = unif_rand()*cumulative_priors[nb_attributes-1]; /* in case priors are not normalized */

    if(nb_attributes == 1) { PRESELECTED_ATTRIBUTES[ith_preselected] = 0; continue;}

// Avoid rounding problems:
    if(x_rand > cumulative_priors[nb_attributes-2]) {
      PRESELECTED_ATTRIBUTES[ith_preselected] = nb_attributes-1;
      continue;
    }
    
    if(x_rand <= cumulative_priors[0]){
      PRESELECTED_ATTRIBUTES[ith_preselected] = 0;
      continue;
    }
    while(lleft != rright -1){
    int mid = (int) ((rright + lleft)/2); /* at least left + 1 */
      if(cumulative_priors[mid] < x_rand){
        lleft = mid;
      }
      else { rright = mid; }
    }
    PRESELECTED_ATTRIBUTES[ith_preselected] = rright;
  }

  /* on construit un noeud et on le place sur le stack */
  index_nodes++;
  prediction[index_nodes]=-1; /* par defaut, c'est pas encore une feuille */
  tested_attribute[index_nodes]=-1; /* pas encore de test */
  left_successor[index_nodes]=-1;
  right_successor[index_nodes]=-1;
  tree=index_nodes;
  
  index_stack_open_nodes++;   /* had been set to -1 before one calls this function? */
  stack_open_nodes[index_stack_open_nodes][0]=tree;
  stack_open_nodes[index_stack_open_nodes][1]=0;
  stack_open_nodes[index_stack_open_nodes][2]=current_learning_set_size-1;
  int current_num_nodes_in_tree = 1;

  /* on lance la boucle de developpement, tant qu'il y a des noeuds ouverts */
  while (index_stack_open_nodes>=0) {
    int node=stack_open_nodes[index_stack_open_nodes][0];
    int start=stack_open_nodes[index_stack_open_nodes][1];
    int end=stack_open_nodes[index_stack_open_nodes][2];
    int nodesize=end-start+1;
    
    /* resume du vecteur */
// wtot = table_score_multiregr[0][2] - (table_score_multiregr[0][1]*table_score_multiregr[0][1])/table_score_multiregr[0][0];
// printf("before summarize_vector: start = %d  end = %d LS_size = %d weighted_variance = %e\n",start,end,current_learning_set_size,wtot);
// printf("before summarize_vector: first row of table_score  = %e  %e  %e\n",table_score_multiregr[0][0],table_score_multiregr[0][1],table_score_multiregr[0][2]);
    summarize_vector(current_learning_set, start, end);
// wtot = table_score_multiregr[0][2] - (table_score_multiregr[0][1]*table_score_multiregr[0][1])/table_score_multiregr[0][0];
// printf("after summarize_vector: first row of table_score  = %e  %e  %e  weighted_var = %e\n",table_score_multiregr[0][0],table_score_multiregr[0][1],table_score_multiregr[0][2],wtot);
    /* pour missing_values */
    node_size[node]=table_score_multiregr[0][0];

    /* condition d'arret: */
    if ((nodesize==1) || (nodesize<min_node_size) || ((MAX_NUM_NODES_IN_TREE != -1) && (current_num_nodes_in_tree >= MAX_NUM_NODES_IN_TREE)) || stop_splitting_criterio()) {  // stop_splitting_criterio is executed only here, checks if variance is too small to bother
      /* c'est une feuille */
      prediction[node]=make_leaf_prediction();
      index_stack_open_nodes--;

    } else {

      /* on recherche le meilleur split */
// printf("before split: v_tot = %e\n",v_tot);
      find_a_split(current_learning_set, start, end);
//        printf("after split: v_tot = %e\n",v_tot);
 //     }

      /* Pas vraiment un test statistique, on vérifie juste si best_threshold_score est >=0 */
//      printf("###### keep_split = %d ######\n",keep_split);
      if ( (keep_split != 0) || not_significant_test()) {
	    /* c'est une feuille */
	    prediction[node]=make_leaf_prediction();
	    index_stack_open_nodes--;
      } else {
	    /* on separe les objets sur base */
	    int left,right,borne;

		borne=separate_ls_vector(best_attribute, best_threshold, current_learning_set,start,end);

		nb_tests++;

                current_num_nodes_in_tree +=2;

		/* on cree deux nouveaux noeuds */
		index_nodes++; left=index_nodes;
		index_nodes++; right=index_nodes;
		prediction[left]=-1; prediction[right]=-1;
		tested_attribute[left]=-1; tested_attribute[right]=-1;
		left_successor[left]=-1; left_successor[right]=-1;
		right_successor[left]=-1; right_successor[right]=-1;

		/* on met a jour le noeud courant */
		threshold[node]=best_threshold;
		tested_attribute[node]=best_attribute;
		left_successor[node]=left-node;
		right_successor[node]=right-node;
	
		/* on place les nouveaux noeuds sur la pile */
		/* pas de best_first, on les met betement sur la pile */
	    stack_open_nodes[index_stack_open_nodes][0]=left;
	    stack_open_nodes[index_stack_open_nodes][1]=start;
	    stack_open_nodes[index_stack_open_nodes][2]=borne-1;
	    index_stack_open_nodes++;
	    stack_open_nodes[index_stack_open_nodes][0]=right;
	    stack_open_nodes[index_stack_open_nodes][1]=borne;
	    stack_open_nodes[index_stack_open_nodes][2]=end;

  
      }     
    }
  }

  /* Loop over all out of bag data and get the associated predictions for accumulation */
  int out_of_bag_data_point, out_of_bag_num=0;
  float *pred, this_average=0.0,this_square=0.0;
  float sum_of_squared_errors=0.0;
 
//  printf("current_learning_set array:\n");
//  for(int jj=0; jj<current_learning_set_size; ++jj) {
//    printf("%d  ",current_learning_set[jj]);
//  }
//  printf("current_learning_set_size = %d nb_obj_in_core_table = %d \n",current_learning_set_size,nb_obj_in_core_table);


// printf("Checking get_tree_prediction \n");
  for (int iobj = (nb_obj_in_core_table-nb_obj_in_validation_table); iobj< nb_obj_in_core_table; ++iobj) {
//    out_of_bag_data_point=1;
//    for(int jj=0; jj<current_learning_set_size; ++jj) {
//      if(current_learning_set[jj] == iobj) out_of_bag_data_point = -1;
//    }
//    if(out_of_bag_data_point == 1) { // out of bag data
      this_average += core_table_y[iobj];
      this_square += core_table_y[iobj]*core_table_y[iobj];
      pred = get_tree_prediction_vector(tree,iobj);
//      ++out_of_bag_num;
      sum_of_squared_errors += (core_table_y[iobj]-pred[0])*(core_table_y[iobj]-pred[0]);
//    }
  }
  this_average = this_average/nb_obj_in_validation_table;
  this_square = this_square/nb_obj_in_validation_table;
  SqE_on_validation[CURRENT_itree] = sum_of_squared_errors/nb_obj_in_validation_table;
  SqE_on_validation[CURRENT_itree] = SqE_on_validation[CURRENT_itree]/(this_square-this_average*this_average); // Normalized to variance of validation data

// printf("finished iobj loop\n");

  return tree;
}

int separate_ls_vector_local(int best_attribute, union threshold_type best_threshold, int *ls_vector, int start, int end) {
  
  while (start!=end) {
    while (start!=end && (check_test(best_attribute,
				     getattval(object_mapping[ls_vector[start]],best_attribute), 
				     best_threshold))) {
      start++;
    }
    while (start!=end && !(check_test(best_attribute,
				      getattval(object_mapping[ls_vector[end]],best_attribute), 
				      best_threshold))) {
      end--;
    }
    if (start!=end) { /* on inverse les deux */
      int temp;
      temp=ls_vector[start];
      ls_vector[start]=ls_vector[end];
      ls_vector[end]=temp;
      start++;
    }
  }
  /* ici, on a start=end, on renvoie la borne */

  if (check_test(best_attribute,getattval(object_mapping[ls_vector[start]],best_attribute), best_threshold))
    return (start+1);
  else
    return start;
}

/* fonction de tri (from the numerical recipes in C) */

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp
#define M_QS 7
#define QUICK_SORT_STACK_SIZE 50
#define VAL(a) getattval(object_mapping[a],att)

/* Here M is the size of subarrays sorted by straight insertion and NSTACK is the required auxiliary storage. */

unsigned long *istack[QUICK_SORT_STACK_SIZE];

void quicksort_ls_vector(int *ls_vector, int start, int end, int att) {  /* Note: att is used through the macro VAL */

  /* Sorts an array arr[1..n] into ascending numerical order using the Quicksort algorithm. n is input; 
     arr is replaced on output by its sorted rearrangement. (extrait de numerical recipes in c) */

  int i,ir=end,j,k,l=start;
  int istack[QUICK_SORT_STACK_SIZE];
  int jstack=-1,o;
  int temp;
  float a;
  
  for (;;) {
    /* Insertion sort when subarray small enough.*/ 
    if (ir-l < M_QS) {
      for (j=l+1;j<=ir;j++) {
	o=ls_vector[j];
	a=VAL(o);
	for (i=j-1;i>=l;i--) {
	  if (VAL(ls_vector[i]) <= a)
	    break;
	  ls_vector[i+1]=ls_vector[i];
	}
	ls_vector[i+1]=o; 
      } 
      if (jstack == -1) 
	break; 
      ir=istack[jstack--];  /* Pop stack and begin a new round of partitioning. */ 
      l=istack[jstack--]; 
    } else { 
      k=(l+ir) >> 1; /* Choose median of left, center, and right elements as partitioning element a. */
                     /* Also rearrange so that a[l]<=a[l+1]<=a[ir]. */

      SWAP(ls_vector[k],ls_vector[l+1]);
      if (VAL(ls_vector[l]) > VAL(ls_vector[ir])) {
	SWAP(ls_vector[l],ls_vector[ir]);
      }
      if (VAL(ls_vector[l+1]) > VAL(ls_vector[ir])) {
	SWAP(ls_vector[l+1],ls_vector[ir]);
      }
      if (VAL(ls_vector[l]) > VAL(ls_vector[l+1])) {
	SWAP(ls_vector[l],ls_vector[l+1]);
      }
      i=l+1; /* Initialize pointers for partitioning. */
      j=ir;
      o=ls_vector[l+1];
      a=VAL(o); /* Partitioning element. */ 
      for (;;) { /* Beginning of innermost loop.*/ 
	do i++; while (VAL(ls_vector[i]) < a); /* Scan up to  nd element > a. */ 
	do j--; while (VAL(ls_vector[j]) > a); /* Scan down to  nd element < a. */ 
	if (j < i) 
	  break; /* Pointers crossed. Partitioning complete. */
	SWAP(ls_vector[i],ls_vector[j]); /* Exchange elements. */ 
      } /* End of innermost loop. */
      ls_vector[l+1]=ls_vector[j]; /* Insert partitioning element. */
      ls_vector[j]=o;
      jstack += 2; /* Push pointers to larger subarray on stack, process smaller subarray immediately. */
      if (jstack > QUICK_SORT_STACK_SIZE) {
	/*Stack too small in quicksort.*/
	return;
      }
      if (ir-i+1 >= j-l) { 
	istack[jstack]=ir; istack[jstack-1]=i; ir=j-1; 
      } else { 
	istack[jstack]=j-1; istack[jstack-1]=l; l=i;
      }
    }
  }
}

/* DISPATCH ON ATTRIBUTE TYPE */
/*----------------------------*/

void find_a_threshold(int att, int *ls_vector, int start, int end) {
  /* simple dispatch sur le type de l'attribut */
  if (NUMERICALp(att)) {
    find_a_threshold_num(att, ls_vector, start, end);
  } else if (SYMBOLICp(att)) {
    find_a_threshold_symb(att, ls_vector, start, end);
  }
}

/* CHECK A TEST */
/* ------------ */

int check_test(int att, float val, union threshold_type threshold) {
  if (NUMERICALp(att)) {
    return (val<threshold.f);
  } else {
    return (BITN(threshold,(int)val));
  }
}

/* ARBRES CLASSIQUES */
/* ----------------- */

/* GENERAL */

/* find_the_best_split */
/* recherche le meilleur split sur tous les attributs */

void find_the_best_split(int *ls_vector, int start, int end) {
  int i;

  best_attribute=-1;
  best_threshold_score=-1.0;
  best_threshold_info=-1.0;

  for(i=0; i<nb_attributes; i++) { /* on considere tous les attributs */

    find_a_threshold(attribute_vector[i], ls_vector, start, end);
    
    if ((current_threshold_score>=0.0) && (current_threshold_score>best_threshold_score)) {
      best_threshold_score=current_threshold_score;
      best_threshold_info=current_threshold_info;
      best_threshold=current_threshold;
      best_attribute=attribute_vector[i];
    }
  }
}

/***************
 * EXTRA-TREES *
 ***************/

float random_split_score_threshold=0.1;
int nb_of_random_tests=1;

void find_a_split_at_random_et(int *ls_vector, int start, int end) {
  int random_att_pos, temp, nb_try=0;
  int remaining_att=nb_attributes;

  best_attribute=-1;
  best_threshold_score=-1.0;
  best_threshold_info=-1.0;
 
  nb_try=0;
  do {
    nb_try++;
    random_att_pos=get_random_integer(remaining_att);
    find_a_threshold(attribute_vector[random_att_pos], ls_vector, start, end);

    if (current_threshold_score>best_threshold_score) {
      best_threshold_score=current_threshold_score;
      best_threshold_info=current_threshold_info;
      best_threshold=current_threshold;
      best_attribute=attribute_vector[random_att_pos];
    }

    remaining_att--;
    if (remaining_att!=0) {
      temp=attribute_vector[random_att_pos];
      attribute_vector[random_att_pos]=attribute_vector[remaining_att];
      attribute_vector[remaining_att]=temp;
    }
    
    if (current_threshold_score<0.0)
      /* l'attribut etait constant, ca ne compte pas */
      nb_try--; 
    
  } while ((remaining_att!=0) && (best_threshold_score<random_split_score_threshold) && (nb_try<nb_of_random_tests));
}

/********************
 * RANDOM FORESTS   *
 ********************/

/* implementation de la variante de Breiman pour les arbres aleatoires
 * on selectionne un certain nombre d'attributs aleatoirement dans l'ensemble
 * des attributs candidats pour lesquels on recherche les meilleurs splits
 */

/* Changes in find_the_best_split_among_k:
 * Added priors when doing the split
 * Removed the part that takes away an attribute (i.e., a gene) once it has
 * already been chosen, this speeds up the code by a factor 2
 */

/* This function sets keep_split to 0 (TRUE) if the splitting should be accepted */
void find_the_best_split_among_k(int *ls_vector, int start, int end) {
  int i;
  int remaining_att=nb_attributes;
  int random_att_pos= 0;
  int pre_att[nb_attributes];
  
  for(int i=0;i<nb_attributes;i++){
    pre_att[i] = -1;
  }
  
  double a, tot;
  
  best_attribute=-1;
  best_threshold_score=-1.0;
  best_threshold_info=-1.0;
  
  int lleft, rright;
  int mid;
  
  for (i=0; (i<rf_k)&&(remaining_att!=0) ; i++) {
    
    tot=0.;
    for(int u=0;u<remaining_att;u++){
      tot+=priors[attribute_vector[u]];
    }
    
    lleft = 0;
    rright = remaining_att;
    
    a = unif_rand()*tot;
    
    if(a<cumulative_priors[0]){
      random_att_pos = lleft;
      
    }else{
      while(lleft != rright -1){
        mid = (int) ((rright + lleft)/2); /* at least left + 1 */
        if(cumulative_priors[mid] <= a){
          lleft = mid;
        }else if(cumulative_priors[mid] >= a){
           rright = mid;
        }
      }
      random_att_pos = rright;
    }
    
    int fflag=0;
    
    for(int j=0;j<remaining_att;j++){
      if(pre_att[j] != random_att_pos){
        fflag++;
      }
    }
    
    if(fflag == remaining_att){
      find_a_threshold(attribute_vector[random_att_pos], ls_vector, start, end);
    }
    
   
    if (current_threshold_score>best_threshold_score) {
      best_threshold_score=current_threshold_score;
      best_threshold_info=current_threshold_info;
      best_threshold=current_threshold;
      best_attribute=attribute_vector[random_att_pos];
    }
    
    /* no more swapping of the attributes */
    remaining_att--;
    
  }
  
  if((best_threshold_score > MIN_FRAC_OF_VARIANCE_BEFORE_SPLIT) && (best_threshold_score*v_tot > MIN_FRAC_OF_TOTAL_WEIGHTED_VARIANCE_IN_SPLIT*total_weighted_variance_at_tree_root)) { 
    keep_split = 0;
    WEIGHTED_VARIANCE_EXPLAINED += best_threshold_score*v_tot;
  }
  else {keep_split = -1; }

}

/* In this variant, the attributes have been previously selected, and so for each 
 * split we get the best attribute in that fixed preselected list, list that is kept unchanged
 * within a tree but varies from tree to tree.
 */

/* Changes in find_the_best_split_among_k_preselected:
 * No priors needed since the attributes have already been selected.
 * Those preselected attributes are specified in the array PRESELECTED_ATTRIBUTES
 */

/* This function sets keep_split to 0 (TRUE) if the splitting should be accepted */
void find_the_best_split_among_k_preselected(int *ls_vector, int start, int end) {
  int ith_preselected, attribute_number;

  best_attribute=-1;
  best_threshold_score=-1.0;
  best_threshold_info=-1.0;

  /* consider all preselected attributes */
  for (ith_preselected=0; ith_preselected<rf_k; ith_preselected++) {
      attribute_number = PRESELECTED_ATTRIBUTES[ith_preselected];
      find_a_threshold(attribute_vector[attribute_number], ls_vector, start, end);
   
    if ((current_threshold_score>=0.0) && (current_threshold_score>best_threshold_score)) {
      best_threshold_score=current_threshold_score;
      best_threshold_info=current_threshold_info;
      best_threshold=current_threshold;
      best_attribute=attribute_vector[attribute_number];
    }
  }
  
  if((best_threshold_score > MIN_FRAC_OF_VARIANCE_BEFORE_SPLIT) && (best_threshold_score*v_tot > MIN_FRAC_OF_TOTAL_WEIGHTED_VARIANCE_IN_SPLIT*total_weighted_variance_at_tree_root)) { 
    keep_split = 0;
    WEIGHTED_VARIANCE_EXPLAINED += best_threshold_score*v_tot;
  }
  else {keep_split = -1; }
}


/******************************************************************************/
/* FONCTIONS UTILITAIRES DIVERSES */

/* valable uniquement si max_val est faible par rapport a rand_max */
/* sinon, il vaut mieux utiliser get_random_long_long */

int get_random_integer(int max_val) {
  /*return (int)floor((double)rand()*max_val*1.0/(RAND_MAX+1.0));*/

	/* change for the R package */
	double R_random_number = unif_rand()*RAND_MAX;
	return (int)floor((double)R_random_number*max_val*1.0/(RAND_MAX+1.0));
}

float get_random_float() {
  /*return (float)((double)rand()*1.0/(RAND_MAX+1.0));*/

    /* change for the R package */
	double R_random_number = unif_rand()*RAND_MAX;
	return (float)((double)R_random_number*1.0/(RAND_MAX+1.0));
}


/* accede a la table de valeurs */

float getattval(int obj, int att) {
  return (float)core_table[att*nb_obj_in_core_table+obj];
}



/* UTILISATION DES ARBRES */

/****************************************
 * propagation d'un objet dans un arbre *
 ****************************************/

/* application d'un arbre a un objet */
/* returns pointer to array whose zeroth element is the zeroth predictor */

float *get_tree_prediction_vector_classical(int tree, int obj) {
  int current_node=tree;
  while (left_successor[current_node]!=-1) {
    if (check_test(tested_attribute[current_node],
		   getattval(obj,tested_attribute[current_node]), 
		   threshold[current_node]))
      current_node+=left_successor[current_node];
    else
      current_node+=right_successor[current_node];
  }
  return prediction_values[prediction[current_node]];
}


/************************************************************************
 * INTERFACE AVEC LE LISP *
 **************************/

/* efface tous les arbres */

void clean_all_trees() {
  /* efface toutes les tables */
  index_nodes=-1;  // keeps track of number of nodes in the growing tree
  index_prediction_values=-1;
  index_stack_open_nodes=-1;
  clean_all_ensemble();
}

/* alloue toutes les tables de donnÈes de maniËre dynamique */

SCORE_TYPE **allocate_table_score_type(int nbl, int nbc) {
  SCORE_TYPE **tableau;
  int i,j;

  tableau=(SCORE_TYPE **)MyMalloc(nbl*sizeof(SCORE_TYPE *));
  if (tableau==NULL)
    return NULL;
  for (i=0;i<nbl;i++) {
    tableau[i]=(SCORE_TYPE *)MyMalloc(nbc*sizeof(SCORE_TYPE));
    if (tableau[i]==NULL) {
      for (j=0; j<i; j++) {
	MyFree((SCORE_TYPE *)tableau[j]);
      }
      return NULL;
    }
  }
  return tableau;
}

void free_table_score_type(SCORE_TYPE **tableau, int nbl) {
  int i;

  if (tableau!=NULL) {
    for (i=0;i<nbl;i++) {
      if (tableau[i]!=NULL)
	MyFree((SCORE_TYPE *)tableau[i]);
    }
    MyFree((SCORE_TYPE **)tableau);
  }
}

float **allocate_table_float(int nbl, int nbc) {
  float **tableau;
  int i,j;

  tableau=(float **)MyMalloc(nbl*sizeof(float *));
  if (tableau==NULL)
    return NULL;
  for (i=0;i<nbl;i++) {
    tableau[i]=(float *)MyMalloc(nbc*sizeof(float));
    if (tableau[i]==NULL) {
      for (j=0; j<i; j++) {
	MyFree((float *)tableau[j]);
      }
      return NULL;
    }
  }
  return tableau;
}

void free_table_float(float **tableau, int nbl) {
  int i;

  if (tableau!=NULL) {
    for (i=0;i<nbl;i++) {
      if (tableau[i]!=NULL)
	MyFree((float *)tableau[i]);
    }
    MyFree((float **)tableau);
  }
}


void free_tree_tables() {
  clean_all_trees();

  if (left_successor!=NULL) {
    MyFree((int *)left_successor);
    left_successor=NULL;
  }
  if (right_successor!=NULL) {
    MyFree((int *)right_successor);
    right_successor=NULL;
  }
  if (tested_attribute!=NULL) {
    MyFree((int *)tested_attribute);
    tested_attribute=NULL;
  }
  if (prediction!=NULL) {
    MyFree((int *)prediction);
    prediction=NULL;
  }
  if (node_weight!=NULL) {
    MyFree((float *)node_weight);
    node_weight=NULL;
  }
  if (node_size!=NULL) {
    MyFree((float *)node_size);
    node_size=NULL;
  }
  if (threshold!=NULL) {
    MyFree((union threshold_type *)threshold);
    threshold=NULL;
  }

  free_table_float(prediction_values,size_current_tree_table_pred);
  prediction_values=NULL;

  free_table_score_type(table_score,3);
  table_score=NULL;

  free_table_score_type(table_score_symb,MAX_NUMBER_OF_SYMBOLIC_VALUES);
  table_score_symb=NULL;

  free_multiregr_table_score();

  size_current_tree_table=0;
  size_current_tree_table_pred=0;

}

int allocate_tree_tables(int nb_of_nodes, int nb_of_leaves, int nb_pred, int tsp) {
  free_tree_tables();

  left_successor=(int *)MyMalloc((size_t)nb_of_nodes*sizeof(int));
  if (left_successor==NULL) {
    return 0;
  }
  right_successor=(int *)MyMalloc((size_t)nb_of_nodes*sizeof(int));
  if (right_successor==NULL) {
    free_tree_tables();
    return 0; 
  }
  tested_attribute=(int *)MyMalloc((size_t)nb_of_nodes*sizeof(int));
  if (tested_attribute==NULL) {
    free_tree_tables();
    return 0;
  }

  node_weight=(float *)MyMalloc((size_t)nb_of_nodes*sizeof(float));
  if (node_weight==NULL) {
    free_tree_tables();
    return 0;
  }
  node_size=(float *)MyMalloc((size_t)nb_of_nodes*sizeof(float));
  if (node_size==NULL) {
    free_tree_tables();    
    return 0;
  }
  threshold=(union threshold_type *)MyMalloc((size_t)nb_of_nodes*sizeof(union threshold_type));
  if (threshold==NULL) {
    free_tree_tables();    
    return 0;
  }

  /* ne sert a rien si multiregr_savepred est a 1 */
  prediction=(int *)MyMalloc((size_t)nb_of_nodes*sizeof(int));
  if (prediction==NULL) {
    free_tree_tables();
    return 0;
  }

  nb_of_predictions=nb_pred;

// DEBUG
   prediction_values=allocate_table_float(nb_of_nodes,1);  /* 2 dimensional table, second argument should be nb_goal_multiregr */
//

  if (nb_pred>0) {

    prediction_values=allocate_table_float(nb_of_leaves,nb_pred);
    if (prediction_values==NULL) {
      free_tree_tables();
      return 0;
    }
    
    /* allocation de la table de score (a ne pas faire si multiregr) */
    if (tsp==1) {
      table_score=allocate_table_score_type(3,nb_pred+1);
      if (table_score==NULL) {
	free_tree_tables();
	return 0;
      }
      table_score_symb=allocate_table_score_type(MAX_NUMBER_OF_SYMBOLIC_VALUES,nb_pred+1);
      if (table_score==NULL) {
	free_tree_tables();
	return 0;
      }
    }

    size_current_tree_table_pred=nb_of_leaves;

  } else
    size_current_tree_table_pred=0;

  size_current_tree_table=nb_of_nodes;

  return 1;

}


/* fonction de test par defaut */

float *(*get_tree_prediction_vector)(int tree, int obj)=get_tree_prediction_vector_classical;

/* pour changer la fonction de test */

void set_test_classical() {
  get_tree_prediction_vector=get_tree_prediction_vector_classical;
}




/*********************/
/* ENSEMBLE D'ARBRES */
/*********************/

/* description de l'ensemble d'arbres */
#define MAX_NUMBER_OF_TREES 10000
int ltrees[MAX_NUMBER_OF_TREES];  /* provides the node number for the tree's root in the forest */ 
float ltrees_weight[MAX_NUMBER_OF_TREES];
// float ltrees_nb_nodes[MAX_NUMBER_OF_TREES];  /* not created as one can access it via get_tree_nb_nodes(ltrees(itree))
float ltrees_LS_size[MAX_NUMBER_OF_TREES];      /* Learning Set size (depends on bagging) for each tree */

int current_nb_of_ensemble_terms=0;
int average_predictions_ltrees=1; /* 1 -> divide pred by the sum of weights 0-> no */

/* DiffÈrentes variantes pour la creation du LS */

/* identite */

float make_ls_vector_identity(int tree) {
  /* on ne fait rien */
  return 1.0;
}

/* bagging */

float make_ls_vector_bagging(int tree) {  /* Argument not used ? */
  int i;

  /* on remet les poids a zero */
  for (i=0; i<global_learning_set_size; i++)
    object_weight[i]=0.0;

  /* on incremente les poids d'objets tires au hasard */
  for (i=0; i<global_learning_set_size; i++) {
    int rn=get_random_integer(global_learning_set_size);
    object_weight[rn]+=1.0;
  }
  
  /* on construit le ls en prenant les objets de poids non nuls */
  /* Thus the LS has a size less or equal to the original dataset and a priori there are no duplications */
  current_learning_set_size=0;
  for (i=0; i<global_learning_set_size; i++) {
    if (object_weight[i]!=0.0) {
      current_learning_set[current_learning_set_size]=i;
      current_learning_set_size++;
    }
  }
//  printf("\n Exiting make_ls_vector_bagging: global_learning_set_size = %d current_learning_set_size = %d \n",global_learning_set_size,current_learning_set_size);

  return 1.0;
}


/*----------------------------*/
/* MART regression (Friedman) */
/*----------------------------*/

/* different variants for the testing of trees */

float current_sum_weight=0.0;

void set_ensemble_method_parameters(int i, int nb_terms, int se, int sewg, float mu) {
  if (i==1) {
    make_ls_vector=make_ls_vector_bagging;
    average_predictions_ltrees=1;
  } else {
    make_ls_vector=make_ls_vector_identity;
    average_predictions_ltrees=1;
  }

  store_ensemble=se;
  save_ensemble_while_growing=sewg;

  number_of_ensemble_terms=nb_terms;  
}

/* pour sauvegarder les ls */
int save_ensemble_ls=0;
int *save_ensemble_ls_vector=NULL;
float *save_ensemble_ls_weight=NULL;
int save_ensemble_ls_size[MAX_NUMBER_OF_TREES];
int save_ensemble_ls_pos=0;
int save_ensemble_ls_nb_ls=0;

void init_save_ensemble_ls(int b) {

  if (save_ensemble_ls && (save_ensemble_ls_vector!=NULL)) {
    MyFree((int *)save_ensemble_ls_vector);
    save_ensemble_ls_vector=NULL;
    MyFree((float *)save_ensemble_ls_weight);
    save_ensemble_ls_weight=NULL;
  }

  save_ensemble_ls=b;

  if (save_ensemble_ls) {
    int s=number_of_ensemble_terms*global_learning_set_size;
    save_ensemble_ls_vector=(int *)MyMalloc((size_t)s*sizeof(int));
    save_ensemble_ls_weight=(float *)MyMalloc((size_t)s*sizeof(float));
    save_ensemble_ls_pos=0;
    save_ensemble_ls_nb_ls=0;
  }
}

/* Construction de l'ensemble d'arbres */

float build_one_tree_ensemble(int *ts_vector, int length_ts_vector) {
  int i, t;
  int sum_complexity=0;
  float current_weight; 
  int nbn=0;
  //FILE *fp;
  float FOREST_AVERAGE_VAR_EXPLAINED;

  /* on vide tout */
  clean_all_trees();
  current_nb_of_ensemble_terms=0;
  current_sum_weight=0.0;

  /* verification de la memoire */
  nbn=(2*global_learning_set_size-1);

  if (!store_ensemble) { /* ensemble non stocke */
    if (size_current_tree_table < nbn) {
      /*memoire trop faible pour construire un arbre*/
      return -1.0;
    }
  } else {/* ensemble stocke en memoire totalement */
    
    if (size_current_tree_table <(number_of_ensemble_terms*nbn)) {
      /*memoire trop faible pour construire l'ensemble d'arbres*/
      return -1.0;
    }
    if (number_of_ensemble_terms>MAX_NUMBER_OF_TREES) {
      return -1.0;
    } 
  }

/*  if (save_ensemble_while_growing) {
  fp=fopen("temp-ensemble-trees.dat", "w");   This crushes the file for all previous targets ! 
   fwrite(&average_predictions_ltrees, sizeof(float), 1, fp);
    fprintf(fp,"%s %d\n","Number of trees = ",number_of_ensemble_terms);
    fclose(fp);
  } */
  
  /* allocation de la matrice de test si necessaire */
  /* de la taille (length_ts_vector*nb_pred) */
  
  /* initialisation de l'ensemble d'apprentissage */
  /* It seems it does this (bagging) for the first tree, then in the loop over trees it makes the next one for each subsequent tree? */
  // make_ls_vector(-1);
  
  
  current_weight = make_ls_vector(-1);  

  FOREST_AVERAGE_VAR_EXPLAINED = 0.0;
  /* boucle de construction */
  for (t=0; t<number_of_ensemble_terms; t++) {
    CURRENT_itree= t;
    int current_tree;  /* actually the node in the forest that stores this tree's root */

    /* si on le demande, on sauve les LS */
    /* pour compute_node_subset_current_ensemble. Idealement, on devrait l'implementer
     * aussi pour le boosting en tenant compte des poids. Ca permettrait de faire des calculs
     * de variable importance plus precise.
     */
    if (save_ensemble_ls) {
      save_ensemble_ls_size[save_ensemble_ls_nb_ls]=current_learning_set_size;
      save_ensemble_ls_nb_ls++;
      for (i=0; i<current_learning_set_size; i++) {
	save_ensemble_ls_vector[save_ensemble_ls_pos]=current_learning_set[i];
	save_ensemble_ls_weight[save_ensemble_ls_pos]=object_weight[current_learning_set[i]];
	save_ensemble_ls_pos++;
      }
    }
    /* Save the LS size */
    ltrees_LS_size[t] = current_learning_set_size;

    /* construction du modele */
    current_tree=build_one_tree();
    
    /*fp=fopen("temp-ensemble-trees.dat", "a");
    fprintf(fp,"itree = %d  fraction_of_weighted_variance_explained = %f  num_nodes = %d LS_size = %d\n",
      t,WEIGHTED_VARIANCE_EXPLAINED/total_weighted_variance_at_tree_root,
      get_tree_nb_nodes(current_tree),current_learning_set_size);
//    ocm_write_one_tree(current_tree, fp);
    fclose(fp);*/
    
    FOREST_AVERAGE_VAR_EXPLAINED += WEIGHTED_VARIANCE_EXPLAINED/total_weighted_variance_at_tree_root;

    sum_complexity+=index_nodes-current_tree+1;  /* index_nodes-current_tree+1 is the number of nodes in the tree, same as get_tree_nb_nodes(current_tree) */
    current_weight=make_ls_vector(current_tree); /* This is for the next tree ! */

    /* stockage du modele si on ne teste pas tout de suite */
    /*if (save_ensemble_while_growing) {
      fwrite(&current_weight,sizeof(float),1,fp);
      write_one_tree(current_tree, fp);
    }*/

    if (store_ensemble) {
      ltrees[t]=current_tree;
      ltrees_weight[t]=current_weight;
      current_nb_of_ensemble_terms++;
    } else
      clean_all_trees();
    
    if (current_weight==0.0) {
      /* on arrete, le modele precedent est trop mauvais */
      t=number_of_ensemble_terms;
      printf("unable to construct tree");
    } else if (current_weight<0) {
      /* on arrete, le dernier modele est parfait. Il devient le seul
	 modele */
      t=number_of_ensemble_terms;
      printf("unable to construct tree");
      if (store_ensemble) {
	current_nb_of_ensemble_terms=1;
	ltrees[0]=ltrees[t];
	ltrees_weight[0]=1.0;
      }
      sum_complexity=index_nodes-current_tree+1;
    }
    
  }
  
  /*if (save_ensemble_while_growing) {
    fclose(fp);
    fp=fopen("temp-ensemble-nb-trees.dat", "wb");
    fwrite(&current_nb_of_ensemble_terms, sizeof(int), 1, fp);
    fclose(fp);
  }*/

  // *FOREST_QUALITY = FOREST_AVERAGE_VAR_EXPLAINED/number_of_ensemble_terms;
  return -1.0;
}


/* pour virer l'ensemble */

void clean_all_ensemble () {
  current_nb_of_ensemble_terms=0;
}


/****************************************
 * SAUVEGARDE ET CHARGEMENT DES MOD»LES *
 ****************************************/

/* calcule le nombre de noeud dans un arbre */

int get_tree_nb_nodes(int tree) {
  int nb_nodes=1;
  int current_node;
  
  index_stack_open_nodes=-1;
  index_stack_open_nodes++;   /* current as well as the last entry of stack_open_nodes, under study */
  stack_open_nodes[index_stack_open_nodes][0]=tree;  /* growing and decreasing stack that contains the node numbers */
  
  while(index_stack_open_nodes>=0) {
    current_node=stack_open_nodes[index_stack_open_nodes][0];
    index_stack_open_nodes--;
    if (left_successor[current_node]!=-1) {
      nb_nodes+=2;
      index_stack_open_nodes++;
      stack_open_nodes[index_stack_open_nodes][0]=current_node+left_successor[current_node];
      index_stack_open_nodes++;
      stack_open_nodes[index_stack_open_nodes][0]=current_node+right_successor[current_node];
    }
  }

  return nb_nodes;

}

/* Ecrit un arbre dans le fichier (version propre a l'algorithme de construction) */

void ocm_write_one_tree(int tree, FILE *fp) {
  int nb_nodes=get_tree_nb_nodes(tree);
  int current_node;

  /* on ecrit le nombre de noeuds */
  // fwrite(&nb_nodes, sizeof(int), 1, fp);
  //fprintf(fp,"nb_nodes = %d \n",nb_nodes);

  /* on ecrit l'info sur les noeuds */
  for (current_node=tree; current_node<tree+nb_nodes; current_node++) {
    int pred;
    /* on ecrit les info sur ce noeud dans le fichier */
    /* on ecrit -1 ou 1 selon la valeur de prediction[current_node] */
    /* deux cas: */
    if (left_successor[current_node]!=-1) {
      pred=-1;
      // fwrite(&pred, sizeof(int), 1, fp);
  //fprintf(fp,"pred = %d  ",pred);
      // fwrite(&tested_attribute[current_node], sizeof(int), 1, fp);
  //fprintf(fp,"tested_attribute = %d  ",tested_attribute[current_node]);
      // fwrite(&threshold[current_node], sizeof(union threshold_type), 1, fp);
  //fprintf(fp,"threshold = %f  ",threshold[current_node].f);
      // fwrite(&left_successor[current_node], sizeof(int), 1, fp);
  //fprintf(fp,"left_successor = %d  ",current_node+left_successor[current_node]);
      /* normalement, on n'a pas besoin de cette valeur 
      // fwrite(&right_successor[current_node], sizeof(int), 1, fp);
  fprintf(fp,"right_successor = %d  \n",current_node+right_successor[current_node]);*/
    } else {
      pred=1;
      // fwrite(&pred, sizeof(int), 1, fp);
  //fprintf(fp,"pred = %d ",pred);
      // on ecrit la valeur de prediction 
      // fwrite(prediction_values[prediction[current_node]], sizeof(float), nb_classes, fp);
  //fprintf(fp,"prediction_entry = %d prediction_val = %f \n",prediction[current_node],prediction_values[prediction[current_node]][0]);
  // Not clear what predition_entry is 
    }
  
  }
      
}


void write_one_tree(int tree, FILE *fp) {
  int nb_nodes=get_tree_nb_nodes(tree);
  int current_node;

  /* on ecrit le nombre de noeuds */
  //fwrite(&nb_nodes, sizeof(int), 1, fp);

  /* on ecrit l'info sur les noeuds */
  for (current_node=tree; current_node<tree+nb_nodes; current_node++) {
    int pred;
    /* on ecrit les info sur ce noeud dans le fichier */
    /* on ecrit 0 ou 1 selon la valeur de prediction[current_node] */
    /* deux cas: */
    if (left_successor[current_node]!=-1) {
      pred=-1;
      fwrite(&pred, sizeof(int), 1, fp);
      fwrite(&tested_attribute[current_node], sizeof(int), 1, fp);
      fwrite(&threshold[current_node], sizeof(union threshold_type), 1, fp);
      fwrite(&left_successor[current_node], sizeof(int), 1, fp);
      /* normalement, on n'a pas besoin de cette valeur */
      fwrite(&right_successor[current_node], sizeof(int), 1, fp);
    } else {
      pred=1;
      fwrite(&pred, sizeof(int), 1, fp);
      /* on ecrit les valeurs de predictions */
      fwrite(prediction_values[prediction[current_node]], sizeof(float), nb_classes, fp);
    }
  }
}







/**************************************************/
/* From here: what used to be in tree-multiregr.c */
/**************************************************/




/*
 * Author: Pierre Geurts (p.geurts@ulg.ac.be)
 * (c) 2002-2010
 *
 * (Please do not distribute without permission of the author)
 *
 */


/******************************/
/* ALGORITHME DE CONSTRUCTION */
/******************************/

int multiregr_savepred=0;

float (*getobjy_multiregr_learn)(int obj, int i); /* permet de separer la table des entrees de la table des sorties */

int allocate_multiregr_table_score(int nb_goal) {

  table_score_multiregr=allocate_table_score_type(3,2*nb_goal+1);
  if (table_score_multiregr==NULL) {
    free_tree_tables();
    return 0;
  }
  
  table_score_symb_multiregr=allocate_table_score_type(MAX_NUMBER_OF_SYMBOLIC_VALUES,2*nb_goal+1);
  if (table_score_symb_multiregr==NULL) {
    free_tree_tables();
    return 0;
  }
  return 1;
}

void free_multiregr_table_score() {
  int i;

  if (table_score_multiregr!=NULL) {
    for (i=0;i<3;i++) {
      if (table_score_multiregr[i]!=NULL)
	MyFree((SCORE_TYPE *)table_score_multiregr[i]);
    }
    MyFree((SCORE_TYPE **)table_score_multiregr);
    table_score_multiregr=NULL;
  }
  if (table_score_symb_multiregr!=NULL) {
    for (i=0;i<MAX_NUMBER_OF_SYMBOLIC_VALUES;i++) {
      if (table_score_symb_multiregr[i]!=NULL)
	MyFree((SCORE_TYPE *)table_score_symb_multiregr[i]);
    }
    MyFree((SCORE_TYPE **)table_score_symb_multiregr);
    table_score_symb_multiregr=NULL;
  }
}

int *goal_multiregr;
int nb_goal_multiregr;  /* Here, it is set to 1 everywhere */

/* calcule les stats sur les sorties */
/* Called just once right after LS is defined */
void summarize_vector_multiregr(int *vector, int start, int end) {
  int i,j;
  SCORE_TYPE w;
  
  
  for (i=0; i<2*nb_goal_multiregr+1; i++) {
    table_score_multiregr[0][i]=0.0;
  }

  for (i=start; i<= end; i++) {
    w=object_weight[vector[i]];
    table_score_multiregr[0][0]+=w;
    for (j=0; j<nb_goal_multiregr; j++) {
      float y=getobjy_multiregr_learn(vector[i],j);
      table_score_multiregr[0][2*j+1]+=w*y;
      table_score_multiregr[0][2*j+2]+=w*y*y;
    }
  }
}

/* fonction de repropagation des objets dans l'arbre */

/* SCORE_TYPE v_tot; */

/* le score est simplement la somme des réductions de variance */
/* Actually, the returned value seems to be reduction_of_weighted_variance/initial_weighted_variance */
/* thus v_tot would be the variance before and "info" would be the reduction_of_variance */

// v_tot is tot_mass*tot_variance before split
// y_tot is y_mass*y_variance where y refers to yes on condition below threshold
// n_tot is n_mass*n_variance where n refers to no on condition below threshold
// v_tot = y_tot + n_tot + y_mass*(mu_y - mu_tot)**2 + n_mass*(mu_n - mu_tot)**2
// This function thus returns the reduction of weighted variance divided by the initial variance

/* For one split, that is the weighted average of (mu_1-mu)^2 and of (mu_2-mu)^2 */
/* where mu is the mean before the split and mu_1 and mu_2 are the means after the split */
/* while the weights are the proportions of data points arising in each part */

SCORE_TYPE compute_multiregr_score_from_table() {
  int i;
  SCORE_TYPE y_tot_var, n_tot_var;

  y_tot_var=0.0;
  n_tot_var=0.0;

  table_score_multiregr[2][0]=table_score_multiregr[0][0]-table_score_multiregr[1][0];  // First column is number of data points = weight
  
  for (i=0; i<nb_goal_multiregr; i++) {
    table_score_multiregr[2][2*i+1]=table_score_multiregr[0][2*i+1]-table_score_multiregr[1][2*i+1];  // sum of all values - sum of all y_values
    table_score_multiregr[2][2*i+2]=table_score_multiregr[0][2*i+2]-table_score_multiregr[1][2*i+2];  // sum of all value**2 - sum of all y_values**2
    
    y_tot_var+=fabs(table_score_multiregr[1][2*i+2]-(table_score_multiregr[1][2*i+1]*table_score_multiregr[1][2*i+1])/table_score_multiregr[1][0]);
    n_tot_var+=fabs(table_score_multiregr[2][2*i+2]-(table_score_multiregr[2][2*i+1]*table_score_multiregr[2][2*i+1])/table_score_multiregr[2][0]);
  }
  
  info=v_tot-(y_tot_var+n_tot_var);
// DEBUG
// printf("In score calculation: num_data_points_at_node = %f  v_tot = %e  y_tot = %e  n_tot = %e  \n",table_score_multiregr[0][0],v_tot,y_tot_var,n_tot_var);

  return (info/v_tot); /* le score est compris entre 0 et 1. La
			  normalisation n'est pas nécessaire */
}

/* stop splitting criterio pour la distance */

/* SCORE_TYPE v_min=0.0; */
/* verifier ce calcul) */

// int stop_splitting_criterio_multiregr_old() {
//   int i;

//  v_tot=0.0;
//  for (i=0; i<nb_goal_multiregr; i++) {
//    v_tot+=table_score_multiregr[0][2*i+2]-(table_score_multiregr[0][2*i+1]*table_score_multiregr[0][2*i+1])/table_score_multiregr[0][0];
//  }
//  return ((v_tot/table_score_multiregr[0][0])<=v_min);
//}

int stop_splitting_criterio_multiregr() {  /* calculates the weighted variance (Sum w_i*y_i**2 - (Sum_i w_i y_i)**2/Sum_i w_i) */
                                           /* denominator is the weighted number of data points (e.g. with bagging, the weights are 1, 2, ...). */
// executed for each splitting node to compute v_tot and check criterion for trying to split */
   int i;
   float ratio;

  v_tot=0.0;
  for (i=0; i<nb_goal_multiregr; i++) {
    v_tot+=table_score_multiregr[0][2*i+2]-(table_score_multiregr[0][2*i+1]*table_score_multiregr[0][2*i+1])/table_score_multiregr[0][0];
  }
// DEBUG
// Check the v_tot calculated here
//  printf("In stop: weighted variance at node before splitting = %e \n",v_tot);
  return ((v_tot/table_score_multiregr[0][0])<=v_min);  // Criterion divides by effective number of data points considered
}


int not_significant_test_multiregr() {
  if (best_threshold_score>=0.0) {
	  return 0;
  } else {
    return 1;
  }
}


/* on n'enregistre rien pour le moment, la prédiction sera calculée au moment du test. Ca evite de retenir un grand vecteur 
 * pour chaque feuille (mais ça augmente les temps de calcul pour le test)
 */

int make_leaf_prediction_multiregr_nosave() {
  return -1;
}

int make_leaf_prediction_multiregr_savepred() {
  int i;
  index_prediction_values++;

// DEBUG: here and below, conclusion is that the memory for prediction_values has not been set up
//  printf("in make_leaf_prediction_multiregr_savepred, nb_goal_multiregr = %d index_prediction_values = %d      ",nb_goal_multiregr,index_prediction_values); 

//

  for (i=0;i<nb_goal_multiregr; i++) {
//    printf(" -> %f  ",(float)(table_score_multiregr[0][2*i+1]/table_score_multiregr[0][0]));
  prediction_values[index_prediction_values][i]=(float)(table_score_multiregr[0][2*i+1]/table_score_multiregr[0][0]);
  }
  return index_prediction_values;
}

/* recherche le meilleur seuil pour la distance */

/* attribut numérique */

void find_the_best_threshold_multiregr(int att, int *ls_vector, int start, int end) {  /* This is the one used for numerical data */
  float old_val, new_val;
  SCORE_TYPE best_score=-1.0, best_info, current_score, w;
  float best_threshold;
  int st=start,i;

  /* initialisation de la table */
  table_score_multiregr[1][0]=0.0;
  for (i=0; i<nb_goal_multiregr; i++) {
    table_score_multiregr[1][2*i+1]=0.0;
    table_score_multiregr[1][2*i+2]=0.0;
  }
  
  /* on trie l'ensemble selon l'attribut */
  sort_ls_vector(ls_vector, start, end, att);

  /* on parcourt toutes les valeurs de seuils possibles */
  old_val=getattval(object_mapping[ls_vector[start]],att);
  
  for(st=start; st<end; st++) {
    w=object_weight[ls_vector[st]];
    table_score_multiregr[1][0]+=w;
    for (i=0; i<nb_goal_multiregr; i++) {
      float y=getobjy_multiregr_learn(ls_vector[st],i);
      table_score_multiregr[1][2*i+1]+=w*y;
      table_score_multiregr[1][2*i+2]+=w*(y*y);
    }

    if ((new_val=getattval(object_mapping[ls_vector[st+1]],att))!=old_val) { /* un nouveau seuil a considerer */
      
      current_score=compute_score_from_table();
      if (current_score>best_score) {
	best_score=current_score;
	best_info=info;
	best_threshold=(old_val+new_val)/2.0;
	if (old_val>=best_threshold) /* problem d'arrondi */
	  best_threshold=new_val;
      }
      old_val=new_val;
    }
  }
  if (best_score>=0.0) {
    current_threshold.f=best_threshold;
    current_threshold_score=best_score;
    current_threshold_info=best_info;
  } else {
    current_threshold_score=-1.0;
  }
}
/* attribut symbolique */

void summarize_symb_att_multiregr(int att, int *vector, int start, int end) {
  int i,j;

  /* set to zero */
  for (i=0; i<NB_VAL(att); i++) {
    table_score_symb_multiregr[i][0]=0.0;
    for (j=0; j<nb_goal_multiregr; j++) {
      table_score_symb_multiregr[i][2*j+1]=0.0;
      table_score_symb_multiregr[i][2*j+2]=0.0;
    }
  }
  
  /* fill the table with frequency */
  for (i=start; i<=end; i++) {
    SCORE_TYPE w=object_weight[vector[i]];
    int v=(int)getattval(object_mapping[vector[i]],att);
    table_score_symb_multiregr[v][0]+=w;
    for (j=0; j<nb_goal_multiregr; j++) {
      float y=getobjy_multiregr_learn(vector[i],j);
      table_score_symb_multiregr[v][2*j+1]+=w*y;
      table_score_symb_multiregr[v][2*j+2]+=w*(y*y);
    }
  }
}

/* ATTENTION version super inefficace et naive */

void find_the_best_threshold_symb_multiregr(int att, int *ls_vector, int start, int end) {
  int i, v, eff_v;
  int nb_val=NB_VAL(att);
  int nb_val_ls=0;
  union threshold_type current_subset, best_subset;
  SCORE_TYPE best_score=-1.0, best_info, current_score;

  /* on precalcule la table avec tous les frequences pour toutes les classes */
  summarize_symb_att_multiregr(att, ls_vector, start, end);

  /* Check that all elements do not belong to the same class */
  for (i=0; i<nb_val; i++) {
    if (table_score_symb_multiregr[i][0]!=0) 
      nb_val_ls++;
  }
  if (nb_val_ls==1) { /* all objects have the same value of this attribute */
    current_threshold_score=-1.0;
    return;
  }

  init_threshold_type(&current_subset);
  add1_threshold_type(&current_subset);

  do {
    /* fill the table score according to the current subset */
    table_score_multiregr[1][0]=0.0;
    for (i=0; i<nb_goal_multiregr; i++) {
      table_score_multiregr[1][2*i+1]=0.0;
      table_score_multiregr[1][2*i+2]=0.0;
    }
      
    eff_v=0;
    for (v=0; v<nb_val; v++) {
      if (table_score_symb_multiregr[v][0]!=0.0) {
	/* check bit eff_v in current_subset */
	if (BITN(current_subset, eff_v)) {
	  table_score_multiregr[1][0]+=table_score_symb_multiregr[v][0];
	  for (i=0; i<nb_goal_multiregr; i++) {
	    table_score_multiregr[1][2*i+1]+=table_score_symb_multiregr[v][2*i+1];
	    table_score_multiregr[1][2*i+2]+=table_score_symb_multiregr[v][2*i+2];
	  }
	}
	eff_v++;
      }
    }
 
    /* compute the score */
    current_score=compute_score_from_table();

    if (current_score>best_score) {
      best_score=current_score;
      best_info=info;
      best_subset=current_subset;
    }
    add1_threshold_type(&current_subset);
  } while (!(BITN(current_subset,(nb_val_ls-1))));
	 
  if (best_score>=0.0) {
    current_threshold_score=best_score;
    current_threshold_info=best_info;
    /* translate current_subset into a proper subset */
    init_threshold_type(&current_threshold);
    eff_v=0;
    for (v=0; v<nb_val; v++) {
      if (table_score_symb_multiregr[v][0]!=0.0) {
	if (BITN(best_subset, eff_v)) {
	  SET_BITN(current_threshold,v);
	}
	eff_v++;
      }
    }    
  } else {
    current_threshold_score=-1.0;
  }
}

/* version extra-trees */

/* attribut numérique */

void find_a_threshold_at_random_multiregr(int att, int *ls_vector, int start, int end) {
  int i, j;
  float min=getattval(object_mapping[ls_vector[start]],att);
  float max=min;
  SCORE_TYPE w;
  
  current_threshold_score=-1.0;

  /* calcule les stats sur l'attribut */
  for (i=start+1; i<=end; i++) {
    float val=getattval(object_mapping[ls_vector[i]],att);
    if (val<min)
      min=val;
    else if (val>max)
      max=val;
  }
  
  if (min==max) { /* toutes les valeurs sont egales */
    return;
  }
  
  /* tirage du seuil (uniformément entre min et max) */
  current_threshold.f=max-(max-min)*get_random_float();
  
  /* calcul du score */
  table_score_multiregr[1][0]=0.0;
  for (i=0; i<nb_goal_multiregr; i++) {
    table_score_multiregr[1][2*i+1]=0.0;
    table_score_multiregr[1][2*i+2]=0.0;
  }

  for (i=start; i<=end; i++) {
    if (getattval(object_mapping[ls_vector[i]],att)<current_threshold.f) {
      w=object_weight[ls_vector[i]];
      table_score_multiregr[1][0]+=w;
      for (j=0; j<nb_goal_multiregr; j++) {
	float y=getobjy_multiregr_learn(ls_vector[i],j);
	table_score_multiregr[1][2*j+1]+=w*y;
	table_score_multiregr[1][2*j+2]+=w*y*y;
      }
    }
  }
  current_threshold_score=compute_score_from_table();
}


/**************************/
/* INTERFACE AVEC LE LISP */
/**************************/

void init_multiregr_trees(int n_min, float vmin, int savepred) {
  min_node_size=n_min;
  v_min=vmin;
  summarize_vector=summarize_vector_multiregr;
  
  multiregr_savepred=savepred;

  if (savepred==1)
    make_leaf_prediction=make_leaf_prediction_multiregr_savepred;
  else
    make_leaf_prediction=make_leaf_prediction_multiregr_nosave;

  // FOR SAVING predictions: we always save!!
  make_leaf_prediction=make_leaf_prediction_multiregr_savepred;

  stop_splitting_criterio=stop_splitting_criterio_multiregr;
  not_significant_test=not_significant_test_multiregr;
  compute_score_from_table=compute_multiregr_score_from_table;
}



/***********************************************/
/* IMPORTANCE DES VARIABLES (version multiple) */
/***********************************************/

/* code separé par rapport aux arbres classiques mais les deux pourraient être combinés sans problème */
/* pour le moment, ts_vector contient des index dans object_mapping (de 0 à N)
 */

/* il faudrait l'allouer de manière dynamique */
int attribute_position[MAX_NUMBER_OF_ATTRIBUTES];

void compute_multiregr_score_from_table_for_varimp(SCORE_TYPE *vi) {
// v_tot is tot_mass*tot_variance before split
// y_tot is y_mass*y_variance where y refers to yes on condition below threshold
// n_tot is n_mass*n_variance where n refers to no on condition below threshold
// v_tot = y_tot + n_tot + y_mass*(mu_y - mu_tot)**2 + n_mass*(mu_n - mu_tot)**2
// This function thus returns the reduction of weighted variance
  int i;
  SCORE_TYPE y_tot_var, n_tot_var;

  table_score_multiregr[2][0]=table_score_multiregr[0][0]-table_score_multiregr[1][0];

  for (i=0; i<nb_goal_multiregr; i++) {
    y_tot_var=0.0;
    n_tot_var=0.0;
    v_tot=0.0;

    v_tot=table_score_multiregr[0][2*i+2]-(table_score_multiregr[0][2*i+1]*table_score_multiregr[0][2*i+1])/table_score_multiregr[0][0];

    table_score_multiregr[2][2*i+1]=table_score_multiregr[0][2*i+1]-table_score_multiregr[1][2*i+1];
    table_score_multiregr[2][2*i+2]=table_score_multiregr[0][2*i+2]-table_score_multiregr[1][2*i+2];
    
    y_tot_var=fabs(table_score_multiregr[1][2*i+2]-(table_score_multiregr[1][2*i+1]*table_score_multiregr[1][2*i+1])/table_score_multiregr[1][0]);
    n_tot_var=fabs(table_score_multiregr[2][2*i+2]-(table_score_multiregr[2][2*i+1]*table_score_multiregr[2][2*i+1])/table_score_multiregr[2][0]);
    vi[i]=v_tot-(y_tot_var+n_tot_var);
  }
}

void get_vi_multiregr_separate(int *ts_vector, int start, int end, int borne, SCORE_TYPE *vi) {
  int i,j;

  /* summarize_vector */
  for (i=0; i<2*nb_goal_multiregr+1; i++) {
    table_score_multiregr[0][i]=0.0;
    table_score_multiregr[1][i]=0.0;
  }

  for (i=start; i<=end; i++) {
    table_score_multiregr[0][0]++;
    for (j=0; j<nb_goal_multiregr; j++) {
      float y=getobjy_multiregr_learn(ts_vector[i],j);
      table_score_multiregr[0][2*j+1]+=y;
      table_score_multiregr[0][2*j+2]+=y*y;
    }
  }
  
  /* calcule v_tot en fait */
  
  if ((start>=borne)||(borne>end)) {
    /* tous les objets sont a droite ou a gauche -> vi=0 
     * (ne peut arriver qu'a cause d'erreur d'arrondi ?)
     */
    for (i=0; i<nb_goal_multiregr; i++) {
      vi[i]=0;
    }
    return;
  }

  /* fill the table_score, il y a sûrement une manière plus efficace que de le faire
   * objet par objet */

  for (i=start; i<borne; i++) {
    table_score_multiregr[1][0]++;
    for (j=0; j<nb_goal_multiregr; j++) {
      float y=getobjy_multiregr_learn(ts_vector[i],j);
      table_score_multiregr[1][2*j+1]+=y;
      table_score_multiregr[1][2*j+2]+=(y*y);
    }
  }

  /* compute the score */
  compute_multiregr_score_from_table_for_varimp(vi);

}

int compute_one_tree_variable_importance_multiregr_separate(int tree, int *ts_vector, int length_ts_vector, float weight,
							     SCORE_TYPE *attribute_importance, int obj) {
  int i, ii;
  SCORE_TYPE *vi=NULL;
  float *Participation_Weights;

  Participation_Weights = (float *)malloc((size_t)(nb_attributes)*sizeof(float));
  for(ii=0;ii<nb_attributes;++ii) Participation_Weights[ii]=0.0;

  vi=(SCORE_TYPE *)malloc((size_t)nb_goal_multiregr*sizeof(SCORE_TYPE));

  if (vi==NULL) {
    /*Impossible d'allouer de la memoire*/
    return 0;
  } else {
  	

	  index_stack_open_nodes++;
	  stack_open_nodes[index_stack_open_nodes][0]=tree;
	  stack_open_nodes[index_stack_open_nodes][1]=0;
	  stack_open_nodes[index_stack_open_nodes][2]=length_ts_vector-1;
  
	  while (index_stack_open_nodes>=0) {
	    int node=stack_open_nodes[index_stack_open_nodes][0];
	    int start=stack_open_nodes[index_stack_open_nodes][1];
	    int end=stack_open_nodes[index_stack_open_nodes][2];
	    int node_size=end-start+1;

	    if ((left_successor[node]==-1)||(node_size==1)) {
	      index_stack_open_nodes--;
	    } else {
	      /* separation */
	      int borne=separate_ls_vector_local(tested_attribute[node], threshold[node], ts_vector, start, end);

	      /* calcul de l'importance (seulement si borne OK) */
	      get_vi_multiregr_separate(ts_vector, start, end, borne, vi);

	      /* mis a jour du vecteur */
	      for (i=0; i<nb_goal_multiregr; i++)
		attribute_importance[i*nb_attributes+attribute_position[tested_attribute[node]]]+=(weight*vi[i]);
                Participation_Weights[tested_attribute[node]] += vi[0];  /* Uses only one goal */
// //               printf("## attribute = %d vi = %e    ",tested_attribute[node],vi[0]);

	      /* left and right successors are put on the stack */
	      index_stack_open_nodes--;
	      if (obj<0) {
		if (start<borne) {
		  index_stack_open_nodes++;
		  stack_open_nodes[index_stack_open_nodes][0]=node+left_successor[node];
		  stack_open_nodes[index_stack_open_nodes][1]=start;
		  stack_open_nodes[index_stack_open_nodes][2]=borne-1;
		}
		if (borne<=end) {
		  index_stack_open_nodes++;
		  stack_open_nodes[index_stack_open_nodes][0]=node+right_successor[node];
		  stack_open_nodes[index_stack_open_nodes][1]=borne;
		  stack_open_nodes[index_stack_open_nodes][2]=end;
		}
	      } else {
		if (check_test(tested_attribute[node],getattval(obj,tested_attribute[node]),threshold[node])) {
		  /* on met le gauche */
		  index_stack_open_nodes++;
		  stack_open_nodes[index_stack_open_nodes][0]=node+left_successor[node];
		  stack_open_nodes[index_stack_open_nodes][1]=start;
		  stack_open_nodes[index_stack_open_nodes][2]=borne-1;
		} else {
		  /* on met le droit */
		  index_stack_open_nodes++;
		  stack_open_nodes[index_stack_open_nodes][0]=node+right_successor[node];
		  stack_open_nodes[index_stack_open_nodes][1]=borne;
		  stack_open_nodes[index_stack_open_nodes][2]=end;
		}
	      }
	    }
	  }
// Get the corresponding Neff = Sum_i w_i**2 / (Sum_i w_i)**2
         float sum=0.0, sum2 = 0.0;
         for(ii=0; ii<nb_attributes; ++ii) {
            sum += Participation_Weights[ii];
            sum2 += Participation_Weights[ii]*Participation_Weights[ii];
         }
         if(sum==0.0) { Neff = 0;}
         else { Neff = sum*sum/sum2; }

	  free((SCORE_TYPE *)vi);
	  free((float *)Participation_Weights);
//	  printf("## Neff = %e\n",Neff);
	  return 1;
  
  }

}

int compute_ltrees_variable_importance_multiregr_separate(SCORE_TYPE *attribute_importance, int obj) {
  /* si obj est à -1 -> calcul classique
   * sinon, on ne parcourt que la branche de l'arbre par laquelle l'objet passe 
   * utile surtout pour un ensemble d'arbre. Permet de retrouver une sorte de branche
   * de tests (floues ou prototype) correspondant à cet objet.
   *
   *  Also for each tree it computes (done in function compute_one_tree_variable_importance_multiregr_separate )
   *  the global variable Neff, that is the effective number of attributes used in the tree (via IPR)
   * (also weighted by the fraction of explained variance at each split)
   * and puts the answer into the array Effective_Num_Attributes_In_Trees having one entry for each tree
   */

  int i,t,j;
  SCORE_TYPE sum_weight=0.0;
  int *ts_vector=current_learning_set;
  int length_ts_vector=global_learning_set_size;
  int flag_continue;

  /* remplit le vecteur */
  for (i=0; i<length_ts_vector; i++)  
    ts_vector[i]=i;  /* The variance used is that for the non bagged data. */

  for (i=0; i<nb_attributes; i++) {
    for (j=0; j<nb_goal_multiregr; j++) {
      attribute_importance[j*nb_attributes+i]=0.0;
    }
    attribute_position[attribute_vector[i]]=i;
  }

  /* boucle sur les arbres */
  if ((current_nb_of_ensemble_terms==0) && (index_nodes>=0)) {
    /* il y a un arbre mais pas d'ensemble. On calcule l'importance de cet arbre uniquement */
	/* flag_continue is meant to replace the exit(0) that used to be in the function compute_one_tree_variable_importance_multiregr_separate*/
    flag_continue = compute_one_tree_variable_importance_multiregr_separate(0, ts_vector, length_ts_vector, 1.0, attribute_importance, obj);
    Effective_Num_Attributes_In_Trees[0] = Neff;  
    printf("######## No ensemble terms\n");
  } else {
    
    for (t=0; t<current_nb_of_ensemble_terms; t++) {
      //ltrees_weight[t] = exp(-SqE_on_validation[t]/full_var/current_nb_of_ensemble_terms);
      //flag_continue = compute_one_tree_variable_importance_multiregr_separate(ltrees[t], ts_vector, length_ts_vector, ltrees_weight[t], attribute_importance, obj);
      //ltrees_weight[t] = exp(-SqE_on_validation[t]/full_var/current_nb_of_ensemble_terms)*1/(1+Neff);
      
      flag_continue = compute_one_tree_variable_importance_multiregr_separate(ltrees[t], ts_vector, length_ts_vector, ltrees_weight[t], attribute_importance, obj);
      
      Effective_Num_Attributes_In_Trees[t] = Neff;    

	  if (flag_continue == 1) {
        sum_weight+=ltrees_weight[t];
	  }
    }
  }

  if (flag_continue == 1) {

     /* normalizing the values */
     if (average_predictions_ltrees==1){
       for (i=0; i<nb_goal_multiregr*nb_attributes; i++)
         attribute_importance[i]/=sum_weight;
     }
  }
  
  return flag_continue; 

}





/*****************************************/
/* From here: what used to be in rtree.c */
/*****************************************/

CORETABLE_TYPE *core_table_y;

float getobjy_multiregr_learn_R(int obj, int att)
{
  return (float)core_table_y[goal_multiregr[att]*nb_obj_in_core_table+object_mapping[obj]];
}

void set_tree_param(int nm, int et, int rf, int rfk) {

  float varmin=0.000000;  /* This parameter sets the minimum variance for splitting a node vs leaving it as a leaf.    */
                     /* Thus setting it to 0 will lead to having at each leaf just one data point (or identical data)  */
                     /* unless some other criterion interrupts the node splitting.                                     */
                     /* Note: the bagging producing the LS does not produce identical data, it just extracts a subset of the data. */

  /* If no need to save the prediction values at the leaves, set savepred to 0, otherwise set to 1 (not working) */
  int savepred = 0; 
  	
  init_multiregr_trees(nm,varmin,savepred);

  set_test_classical();
  
  if (et==1) {
    find_a_threshold_num=find_a_threshold_at_random_multiregr;
    find_a_threshold_symb=find_the_best_threshold_symb_multiregr;
  } else {
    find_a_threshold_num=find_the_best_threshold_multiregr;
    find_a_threshold_symb=find_the_best_threshold_symb_multiregr;
  }

  if (et==1) {
    find_a_split=find_a_split_at_random_et;
    nb_of_random_tests=rfk;
    random_split_score_threshold=10.0;
    rf_k=rfk;
  } else if (rf==1) {
    find_a_split=find_the_best_split_among_k;
    find_a_split=find_the_best_split_among_k_preselected;
    nb_of_random_tests=rfk;
    random_split_score_threshold=10.0;
    rf_k=rfk;
  } else {
    find_a_split=find_the_best_split;
    random_split_score_threshold=10.0;
  }
}

int mart=0;

void set_ensemble_param(int nbterms, int bs) {

  int method;

  if (bs==1)
    method=1;
  else
    method=0;

  set_ensemble_method_parameters(method,nbterms,1,0,1.0);  /* arguments 3: save_ensemble, argument 4: save_ensemble_while_growing  */

}


// build a tree and return the variable importances

// Input arguments:
//  - nbobj_total: number of objects (including validation)
//  - nbobj_validation: number of objects to keep at the end of the dataset for validation
//  - nbatt: number of attributes (loci, pairs of loci...)
//  - X: table with the input data (rows = objects, columns = attributes)
//  - Y: table with the output data
//  - nm: n_min: minimum number of objects to split a node
//  - et: 1 -> extra trees randomisation
//  - rf: 1 -> random forests randomisation
//  - rfk: parameter mtry of random forests or K in extra-trees
//  - nbterms: number of trees in the ensemble
//  - bs: 1 -> bootstrap sampling
//  - fsrf: 0 -> variable importances based on information 1 -> variable importances based on OOB and permutations
//  - P: array of priors on the interactions for the (single) target gene
// Output arguments (memory must have been generated within the R and associated pointers passed):
//  - variable importance values *vimp 
//  - pForest_Prediction_SS_on_validation: squared error on validation data of Forest predictions (sum over iobj of
//   squared error of the Forest prediction, prediction obtained by averaging across trees)
//  - n_eff: number of effective genes in each tree
//  - square_err_validation: (normalised by variance) square error of each tree calculated on the validation data
//  - ls_size: size of the learning set for each tree

void BuildTreeEns(int *nbobj_total, int *nbobj_validation, int *nbatt, CORETABLE_TYPE *X, CORETABLE_TYPE *Y, int *nm, int *et, int *rf, int *rfk, int *nbterms,
                  int *bs, int *fsrf, SCORE_TYPE *vimp, float *P, float *pForest_Prediction_SS_on_validation, float *N_eff, float *square_err_validation, 
                  float *ls_size)
{
  /* Necessary for using unif_rand(), it gets the state of the RNG initialized by R */
  GetRNGstate();
  // printf("first random number = %f \n",get_random_float());
	
  int i, length_ls_vector, maxnbnodes, flag_continue;

  nb_attributes=*nbatt;
  nb_obj_in_core_table=*nbobj_total;
  nb_obj_in_validation_table=*nbobj_validation;

  /* set the core_table and the accessors */
  core_table=X;

  /* set the attributes */
  attribute_descriptors=(int *)MyMalloc((size_t)nb_attributes*sizeof(int));
  length_attribute_descriptors=nb_attributes;
  
  cumulative_priors = (float *)MyMalloc((size_t)nb_attributes*sizeof(float));

  PRESELECTED_ATTRIBUTES = (float *)MyMalloc((size_t)(*rfk)*sizeof(float));

// printf("success allocating PRESELECTED_ATTRIBUTES, rfk=%d \n",*rfk); 
   
  /* get the priors and the cumulative */
  priors = P;
  
  float SUM = 0.;

  for(int i=0; i<((int) nb_attributes);i++){
    SUM+=priors[i];
    cumulative_priors[i]=SUM;
  }
  
  for (i=0;i<nb_attributes;i++)
    /* All attributes are assumed to be numerical */
    attribute_descriptors[i]=0;   

  attribute_vector=(int *)MyMalloc((size_t)nb_attributes*sizeof(int));
  for (i=0; i<nb_attributes;i++)
    attribute_vector[i]=i; 

  /* Build the LS and the weight vector (init_learning_set) */

  length_ls_vector=(nb_obj_in_core_table - nb_obj_in_validation_table);  /* don't use validation data */
  global_learning_set_size=length_ls_vector;
  current_learning_set_size=length_ls_vector;
  object_mapping=(int *)MyMalloc((size_t)length_ls_vector*sizeof(int));
  for (i=0; i<length_ls_vector; i++)
    object_mapping[i]=i;

  current_learning_set=(int *)MyMalloc((size_t)length_ls_vector*sizeof(int));
  for (i=0; i<length_ls_vector; i++)
    current_learning_set[i]=i;

  object_weight=(SCORE_TYPE *)MyMalloc((size_t)length_ls_vector*sizeof(SCORE_TYPE));
  for (i=0;i<length_ls_vector;i++)
    object_weight[i]=1.0;
  
  /* Initialise regression problem */
  nb_goal_multiregr=1;
  core_table_y=Y;
  goal_multiregr=(int *)MyMalloc((size_t)nb_goal_multiregr*sizeof(int));
  for (i=0;i<nb_goal_multiregr; i++)
    goal_multiregr[i]=i;
  getobjy_multiregr_learn=getobjy_multiregr_learn_R;
  
  /* Setting tree parameters */
  set_tree_param(*nm, *et, *rf, *rfk);
  set_ensemble_param(*nbterms, *bs);
  init_save_ensemble_ls(1);

  /* Allocate memory for trees */
  maxnbnodes=number_of_ensemble_terms*(2*length_ls_vector-1);
  allocate_tree_tables(maxnbnodes,ceil((maxnbnodes+number_of_ensemble_terms)/2),0,0);  /* 3rd argument = 0 -> no saving of prediction_values */ 
  allocate_multiregr_table_score(nb_goal_multiregr);
  clean_all_trees();


  Effective_Num_Attributes_In_Trees = N_eff;
  SqE_on_validation = square_err_validation;
  LS_size = ls_size;

  pForest_Prediction_SS_on_validation_to_Return = pForest_Prediction_SS_on_validation;  /* Sets the pointer so one can return the value of that quantity */
  
  //for(int ii=0;ii<number_of_ensemble_terms;++ii) Effective_Num_Attributes_In_Trees[ii]=0.0;

  /* Learn an ensemble of trees */
  
  build_one_tree_ensemble(NULL, 0);

  /* Re-sort the variables */
  for (i=0; i<nb_attributes;i++)
    attribute_vector[i]=i;
  
  /* Compute variable importances */
  flag_continue = compute_ltrees_variable_importance_multiregr_separate(vimp,-1); 
  
  for (int itree = 0; itree<*nbterms; ++itree) {
    LS_size[itree] = ltrees_LS_size[itree];

  }

// printf("calculation of pForest_Prediction_SS_on_validation_to_Return\n");
  /* Get the Forest prediction on validation data (using average over trees) */
  /* and on the fly, get the square error of that Forest prediction, per validating iobj */
  *pForest_Prediction_SS_on_validation_to_Return = 0.0;
  int valid_counts=0;
  double average_pred_y = 0.0;
  double average_y=0.0, average_y_squared=0.0;
  for(int iobj = (nb_obj_in_core_table-nb_obj_in_validation_table); iobj< nb_obj_in_core_table; ++iobj) {
    average_y += core_table_y[iobj];
    average_y_squared += core_table_y[iobj]*core_table_y[iobj];
    for(int itree = 0; itree<*nbterms; ++itree) {
      average_pred_y += get_tree_prediction_vector(ltrees[itree], iobj)[0];
    }
    average_pred_y = average_pred_y/(*nbterms);
    *pForest_Prediction_SS_on_validation_to_Return += (average_pred_y - core_table_y[iobj])*(average_pred_y - core_table_y[iobj]);
  }
  if(nb_obj_in_validation_table > 0) {
    *pForest_Prediction_SS_on_validation_to_Return = *pForest_Prediction_SS_on_validation_to_Return/nb_obj_in_validation_table;
    average_y = average_y/nb_obj_in_validation_table;
    average_y_squared = average_y_squared/nb_obj_in_validation_table;
    *pForest_Prediction_SS_on_validation_to_Return = *pForest_Prediction_SS_on_validation_to_Return/(average_y_squared - average_y*average_y);
  }

  if (flag_continue == 1) {

    /* Free memory */
    MyFree(attribute_descriptors);
    MyFree(attribute_vector);
    MyFree(goal_multiregr);
    MyFree(object_mapping);
    MyFree(current_learning_set);
    MyFree(cumulative_priors);
    MyFree(PRESELECTED_ATTRIBUTES);

    MyFree(object_weight);
    free_tree_tables();
    MyFree((int *)save_ensemble_ls_vector);
    save_ensemble_ls_vector=NULL;
    MyFree((float *)save_ensemble_ls_weight);
    save_ensemble_ls_weight=NULL;
	
  }
  
  /* Necessary for using unif_rand(), it saves the state of the RNG so it is known by R */
  PutRNGstate();

}




/***************************/
/* Register native routine */
/***************************/

#include <R_ext/Rdynload.h>

static R_CMethodDef cMethods[] = {
	{"BuildTreeEns",  (DL_FUNC) &BuildTreeEns, 18},
	{NULL}
};


void R_init_PdynGENIE3(DllInfo *info) {
	/* Register the .C routine. No .Call(), .Fortran() or .External() routines,so pass those arrays as NULL.*/
	R_registerRoutines(info,cMethods, NULL, NULL, NULL);
}
