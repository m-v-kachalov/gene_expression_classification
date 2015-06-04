#include<stdio.h>
#include<stdlib.h>
#include<math.h>


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


double welch_variance(float *expression_vals, int num_class_0, int num_samples){
	float temp = 0.0;	float var_sum_0 = 0.0;	float var_sum_1 = 0.0;
	for (int i = 0; i < num_class_0; i++) { temp += expression_vals[i]; }
	//printf("temp %f\n", temp);
	float class_0_avg = temp/num_class_0;
	temp = 0.0;
	for (int i = num_class_0; i < num_samples; i++) { temp += expression_vals[i]; }
//	printf("temp %f\n", temp);
	float class_1_avg = temp/(num_samples - num_class_0);
// calculate variance for 
	for (int i = 0; i < num_class_0; i++) {
		var_sum_0 += pow(fabs(class_0_avg - expression_vals[i]),2); }

	for (int i = num_class_0; i < num_samples; i++) {
		var_sum_1 += pow(fabs(class_1_avg - expression_vals[i]),2); }

	float variance_0 = var_sum_0/(num_class_0 - 1);
	float variance_1 = var_sum_1/(num_samples - num_class_0 - 1);
	//printf(" avg 0 %f    avg 1 %f     var 0 %f     var 1 %f \n", class_0_avg, class_1_avg, variance_0, variance_1);
	return ( fabs(class_0_avg - class_1_avg) / 
		 sqrt( variance_0/num_class_0 + variance_1/(num_samples - num_class_0)) );
}

int k_neighbours(unsigned k, sample query, sample * train_samples,
			unsigned num_samples, unsigned num_genes){
	float distances[num_samples];
	for (unsigned i = 0; i < num_samples; i++){
		distances[i] = euclidean(query, train_samples[i], num_genes);
	}
	int k_nghbrs[k];
	for (unsigned i = 0; i < k; i++) {
 		float cur_min = 1e+37; int temp = 0;
		for (unsigned j = 0; j < num_samples; j++){
			if ( distances[j] < cur_min && distances[j] > 0.0){
				cur_min = distances[j]; temp = j; }
		}
		k_nghbrs[i] = temp;
		distances[temp] = -100.0;			
	}
	unsigned cntr = 0;
	for (unsigned i = 0; i < k; i++){
		if (train_samples[k_nghbrs[i]].class == 0) { cntr++; }
	}
	if (cntr >= ( k/2 + 1 )) { return 0;} else {return 1;}
}

int welch_k_neighbours(unsigned k, sample query, sample * train_samples,
		 unsigned num_samples, unsigned welch, rank *welch_rank){
	float distances[num_samples];
	for (unsigned i = 0; i < num_samples; i++){
		distances[i] = welch_euclidean(query, train_samples[i], 
						welch, welch_rank);
	}
	int k_nghbrs[k];
	for (unsigned i = 0; i < k; i++) {
 		float cur_min = 1e+37; int temp = 0;
		for (unsigned j = 0; j < num_samples; j++){
			if ( distances[j] < cur_min && distances[j] > 0.0){
				cur_min = distances[j]; temp = j; }
		}
		k_nghbrs[i] = temp;
		distances[temp] = -100.0;			
	}
	unsigned cntr = 0;
	for (unsigned i = 0; i < k; i++){
		if (train_samples[k_nghbrs[i]].class == 0) { cntr++; }
	}
	if (cntr >= ( k/2 )) { return 0;} else {return 1;}
}

// first argument is training data, 

int main(int argc, char ** argv){
	FILE * fp = fopen(argv[1], "r"); FILE * q_fp = fopen(argv[2], "r");
	int num_genes, num_class_0, num_class_1, num_samples;
	fscanf(fp,"%i\t%i\t%i\n", &num_genes, &num_class_0, &num_class_1);
	num_samples = num_class_0 + num_class_1;
	rank *welch_rank = malloc(sizeof(rank) * num_genes);
	printf("%d\t%d\t%d\n", num_genes, num_class_0, num_class_1);
	sample *t_samples = malloc( sizeof(sample) * num_samples); 
	for (int i = 0; i < num_samples ; i++){
		t_samples[i].genes = malloc(sizeof(float) * num_genes);
		t_samples[i].class = ( i < num_class_0 ? 0 : 1);
	}

// credit to stackoverflow.com/questions/5357184 for parsing tab-delimited, numerical data
	char buffer[5000]; float cur_expression_vals[num_samples];
	unsigned cur_gene = 0;
	while (fgets(buffer, sizeof(buffer), fp)){
    		size_t i = 0;
    		char *next = buffer;
    		while(*next && *next != '\n') { 
			cur_expression_vals[i] = strtof(next, &next);
			//printf("cur expression_vals %.2f\n", cur_expression_vals[i]); 
			t_samples[i].genes[cur_gene] = cur_expression_vals[i++];
		}
		if ( (int)i != num_samples ) { printf(" error at gene %i\n", (int) cur_gene); exit(2); }
//		printf("cur_gene %i", (int) cur_gene);
		welch_rank[cur_gene] = malloc(sizeof(rank));	
		welch_rank[cur_gene]->gene = cur_gene;
		double test = welch_variance(cur_expression_vals, num_class_0, num_samples);
	//	printf("welch %f\n", test);
		welch_rank[cur_gene]->t_score = (float)test;
		cur_gene++;
	}
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
	fscanf(q_fp,"%i\t%i\t%i\n", &num_genes, &num_class_0, &num_class_1);
	int num_queries = num_class_0 + num_class_1;
	printf("%d\t%d\t%d\n", num_genes, num_class_0, num_class_1);
	sample *q_samples = malloc( sizeof(sample) * num_queries); 
	for (int i = 0; i < num_queries ; i++){
		q_samples[i].genes = malloc(sizeof(float) * num_genes);
		q_samples[i].class = -1;
	}
	cur_expression_vals[num_queries];
	cur_gene = 0;
	while (fgets(buffer, sizeof(buffer), q_fp)){
    		size_t i = 0;
    		char *next = buffer;
    		while(*next && *next != '\n') { 
			cur_expression_vals[i] = strtof(next, &next);
			q_samples[i].genes[cur_gene] = cur_expression_vals[i];
			i += 1;
		}
		if ((int)i != num_queries ) { printf("error\n"); exit(1); }
		cur_gene++;
	}
	unsigned k = strtoul(argv[3], NULL, 10);
	unsigned welch = 0;
	if (argc > 4) { welch = strtoul(argv[4], NULL, 10); }
	if (welch != 0) { 
	   for (int i = 0; i < num_queries; i++){
		int result = welch_k_neighbours(k, q_samples[i], t_samples, 
					   num_samples, welch, welch_rank);
		printf("welch query %i result %i\n", i, result);
		}
	}
	else {
	   for (int i = 0; i < num_queries; i++){
		int result = k_neighbours(k, q_samples[i], t_samples, 
					   num_samples, num_genes);
		printf("query %i result %i\n", i, result);
	}	
	}	
}


