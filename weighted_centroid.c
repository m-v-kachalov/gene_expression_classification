#include <stdio.h>
#include <stdlib.h>
#include <math.h>


typedef struct sample { float * genes; int class; } sample;

typedef struct rank { double t_score; unsigned gene; } *rank;


float euclidean(sample a, sample b, unsigned num_genes){
	float run_sum = 0.0;
	for (unsigned i = 0; i < num_genes; i++){
		run_sum = pow( fabs(a.genes[i] - b.genes[i]), 2);
	}
	return sqrt(run_sum);
}

float welch_euclidean(sample a, sample b, unsigned welch, rank *welch_rank){
	float run_sum = 0.0;
	for (unsigned i = 0; i < welch; i++){
		rank cur_rank = welch_rank[i];
		run_sum = pow( fabs(a.genes[cur_rank->gene] 
							- b.genes[cur_rank->gene]), 2);
	}
	return sqrt(run_sum);
}


double welch_variance(float *vals, int num_class_0, int num_samples){
	float temp = 0.0;	float var_sum_0 = 0.0;	float var_sum_1 = 0.0;
	for (int i = 0; i < num_class_0; i++) { temp += vals[i]; }
	float class_0_avg = temp/num_class_0;
	temp = 0.0;
	for (int i = num_class_0; i < num_samples; i++) { temp += vals[i]; }
	float class_1_avg = temp/(num_samples - num_class_0);
	for (int i = 0; i < num_class_0; i++) {
		var_sum_0 += pow(fabs(class_0_avg - vals[i]),2); }
	for (int i = num_class_0; i < num_samples; i++) {
		var_sum_1 += pow(fabs(class_1_avg - vals[i]),2); }
	float variance_0 = var_sum_0/(num_class_0 - 1);
	float variance_1 = var_sum_1/(num_samples - num_class_0 - 1);
	return ( fabs(class_0_avg - class_1_avg) / 
		 sqrt( variance_0/num_class_0 + variance_1/(num_samples - num_class_0)) );
}

int weighted_centroid ( sample query, sample *t_samples, unsigned num_samples, float limit,
  unsigned num_class_0, rank *welch_rank, double avg_tscore, unsigned num_genes ){
	double s_class0 = 0.0;	double s_class1 = 0.0;
	for (unsigned x = 0; x < num_genes; x++){
		float welch_scale = welch_rank[x]->t_score;
		if (welch_scale/avg_tscore < limit) { break; }
		unsigned i = welch_rank[x]->gene;
		double cntr_0 = 0.0; double cntr_1 = 0.0;
		double avg_class0 = 0.0; double avg_class1 = 0.0;
		for (unsigned j = 0; j < num_class_0; j++){
			avg_class0 += t_samples[j].genes[i];
			cntr_0 += fabs(query.genes[i] - t_samples[j].genes[i]);
		}
		for (unsigned j = num_class_0; j < num_samples; j++){
			avg_class1 += t_samples[j].genes[i];
			cntr_1 += fabs(query.genes[i] - t_samples[j].genes[i]);
		}
		cntr_0 = cntr_0/num_class_0;
		cntr_1 = cntr_1/(num_samples - num_class_0);
		float total_dist = fabs(cntr_0) + fabs(cntr_1);
		if (cntr_0 < cntr_1) { float temp_class0 = 
						welch_scale * fabs(cntr_1 - cntr_0)/total_dist;
				     	s_class0 += temp_class0; s_class1 -= temp_class0;	
				     }				
		else 		     { float temp_class1 = 
						welch_scale * fabs(cntr_0 - cntr_1)/total_dist; 
					s_class1 += temp_class1;   s_class0 -= temp_class1;
				     }
	}
	if ( s_class0 > s_class1 ) { return 0; } else { return 1; }
} 

int main(int argc, char ** argv){
	FILE * t_fp = fopen(argv[1],"r");  FILE * q_fp = fopen(argv[2],"r"); 
	float scale_limit = 0.0;
	if (argc > 3) { scale_limit = strtof(argv[3],NULL); }
	int num_genes, num_class_0, num_class_1, num_samples;
	fscanf(t_fp,"%i\t%i\t%i\n", &num_genes, &num_class_0, &num_class_1);
	num_samples = num_class_0 + num_class_1;
	rank *welch_rank = malloc(sizeof(rank) * num_genes);
	printf("%d\t%d\t%d\n", num_genes, num_class_0, num_class_1);
	sample *t_samples = malloc( sizeof(sample) * num_samples); 
	for (int i = 0; i < num_samples ; i++){
		t_samples[i].genes = malloc(sizeof(float) * num_genes);
		t_samples[i].class = ( i < num_class_0 ? 0 : 1);
	}

// credit to stackoverflow.com/questions/5357184 for parsing tab-delimited, numerical data
	char buffer[5000]; float cur_vals[num_samples];
	unsigned cur_gene = 0;
	double avg_tscore = 0.0; double max_tscore = 0.0;
	while (fgets(buffer, sizeof(buffer), t_fp)){
    		size_t i = 0;
    		char *next = buffer;
    		while(*next && *next != '\n') { 
			cur_vals[i] = strtof(next, &next);
			t_samples[i].genes[cur_gene] = cur_vals[i++];
		}
		if ( (int)i != num_samples ) { printf(" error at gene %i\n", (int) cur_gene); exit(2); }
		welch_rank[cur_gene] = malloc(sizeof(rank));	
		welch_rank[cur_gene]->gene = cur_gene;
		double test = welch_variance(cur_vals, num_class_0, num_samples);
		if ( test > max_tscore ) { max_tscore = test; }
		welch_rank[cur_gene]->t_score = (float)test;
		avg_tscore += test;
		cur_gene++;
	}
	avg_tscore = avg_tscore/num_genes;
	printf("done reading in\n");
	printf("max tscore : %.2f\n", max_tscore);
	printf("avg tscore : %.2f\n",avg_tscore);
	//exit(1);
	for (int i = 1; i < num_genes; i++){
		int cur = i; rank temp;
		while ( cur > 0 && welch_rank[cur]->t_score > welch_rank[cur - 1]->t_score ) {
			temp = 	welch_rank[cur];
			welch_rank[cur] = welch_rank[cur-1];
			welch_rank[cur-1] = temp;
			cur--;
		}		
	}
	printf("done welch sort\n");
	int num_0_samples = num_class_0;
	fscanf(q_fp,"%i\t%i\t%i\n", &num_genes, &num_class_0, &num_class_1);
	int num_queries = num_class_0 + num_class_1;
	printf("%d\t%d\t%d\n", num_genes, num_class_0, num_class_1);
	sample *q_samples = malloc( sizeof(sample) * num_queries); 
	for (int i = 0; i < num_queries ; i++){
		q_samples[i].genes = malloc(sizeof(float) * num_genes);
		q_samples[i].class = -1;
	}
	cur_vals[num_queries];
	cur_gene = 0;
	while (fgets(buffer, sizeof(buffer), q_fp)){
    		size_t i = 0;
    		char *next = buffer;
    		while(*next && *next != '\n') { 
			cur_vals[i] = strtof(next, &next);
			q_samples[i].genes[cur_gene] = cur_vals[i++];
		}
		cur_gene++;
	}
	for (int i = 0; i < num_queries; i++){
		int result = weighted_centroid(q_samples[i], t_samples, num_samples,
		        scale_limit, num_0_samples ,welch_rank, avg_tscore, num_genes);
		printf("query %i result %i\n", i, result);
	}


}





















