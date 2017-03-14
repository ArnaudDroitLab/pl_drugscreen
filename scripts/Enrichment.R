#setwd("C:/Users/Eric/Desktop/Screen drogue/")

all_data = read.table("input/All.txt", sep="\t", quote="\"", header=TRUE, comment.char = "")

# Clean up columns
all_data$Class <- gsub("\n", "", all_data$Class)
all_data$Action <- gsub("\n", "", all_data$Action)
all_data$Percent.Bloating = as.numeric(as.character(all_data$Percent.Bloating))


# Define hit classes
hit_classes = list(
    # A) 0% Bloating - 100% Survie
    a = (all_data$Percent.Bloating == 0) & (all_data$Percent.survival == 1),
    
    # B) 0% Bloating - 75-100% Survie
    b = (all_data$Percent.Bloating == 0) & (all_data$Percent.survival >= 0.75),
    
    # C) >60 oeufs - 100% Survie
    c = (all_data$Number.of.eggs >= 60) & (all_data$Percent.survival == 1),
    
    # D) >60 oeufs - 75-100% Survie
    d = (all_data$Number.of.eggs >= 60) & (all_data$Percent.survival >= 0.75),
    
    # E) 0% Tumeur
    #e = (all_data$Tumor.percent.corrected == 0),
    
    # F) 100% Bloating + 100% Survie
    f = (all_data$Percent.Bloating == 1) & (all_data$Percent.survival == 1),
    
    # G) 100% Bloat + 75-100% Survie
    g = (all_data$Percent.Bloating == 1) & (all_data$Percent.survival >= 0.75))
    
    # H) 0% Survie
    #h = all_data$Percent.survival == 0)

# Perform enrichment of terms present more than 10 times.
drug_classes_to_test = names(table(all_data$Class)[table(all_data$Class) >= 10])

results = data.frame()
for(i in 1:length(drug_classes_to_test)) {
    drug_class = functional_classes_to_test[i]

    drugs_in_class = sum(all_data$Class==drug_class, na.rm=TRUE)
    total_drugs = nrow(all_data)

    for(hit_class in names(hit_classes)) {
        n_hits = sum(hit_classes[[hit_class]], na.rm=TRUE)
        hits_in_drug_class = sum(all_data$Class[hit_classes[[hit_class]]] == drug_class, na.rm=TRUE)
        expected_hits = n_hits * (drugs_in_class / total_drugs)
        pval = phyper(hits_in_drug_class, 
                      drugs_in_class,
                      total_drugs - drugs_in_class,
                      n_hits,
                      lower.tail=FALSE)
        
        
        results = rbind(results, 
                        data.frame(Drug_Class=drug_class,
                                   Hit_class=hit_class,
                                   Drugs_In_Class=drugs_in_class,
                                   Number_of_hits=n_hits,
                                   Hits_in_drug_class=hits_in_drug_class,
                                   Expected=expected_hits,
                                   Enrichment=hits_in_drug_class / expected_hits,
                                   PVal=pval))
    }
}