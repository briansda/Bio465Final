library(MASS)

output_vector_of_vectors
control_subtype_vector
identify control GSM's place in a vector
identify subtype GSM's place in a vector

comparisons<-control_GSM_vector and subtype_GSM_vector

for each gene in file:
	identify control vs. subtype
	for each comparison in comparisons:
		output<-wilcox.test(control, subtype)
		control_subtype_vector<-output.p-value

output vector_of_vectors
