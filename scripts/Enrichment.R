#setwd("C:/Dev/Projects/pl_screen")

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
    e = (all_data$Tumor.presence == 0 & !is.na(all_data$Tumor.presence)),
    
    # F) 100% Bloating + 100% Survie
    f = (all_data$Percent.Bloating == 1) & (all_data$Percent.survival == 1),
    
    # G) 100% Bloat + 75-100% Survie
    g = (all_data$Percent.Bloating == 1) & (all_data$Percent.survival >= 0.75),
    
    # H) 0% Survie
    h = all_data$Percent.survival == 0)


# Given a data-set, a set of drug classes and sets of hits,
# perform category enrichment on all drug classes ~ hit class
# combinations.
test_drug_classes <- function(all_data, drug_classes_to_test, hit_classes) {
    results = data.frame()
    # Loop over all drug classes
    for(i in 1:length(drug_classes_to_test)) {
        # Get the class name and the number of drugs in it.
        drug_class = drug_classes_to_test[i]
        drugs_in_class = sum(all_data$Class==drug_class, na.rm=TRUE)
        total_drugs = nrow(all_data)
    
        # Loop over all hit classes.
        for(hit_class in names(hit_classes)) {
            # Calculate metrics.
            n_hits = sum(hit_classes[[hit_class]], na.rm=TRUE)
            hits_in_drug_class = sum(all_data$Class[hit_classes[[hit_class]]] == drug_class, na.rm=TRUE)
            expected_hits = n_hits * (drugs_in_class / total_drugs)
            pval = phyper(hits_in_drug_class, 
                          drugs_in_class,
                          total_drugs - drugs_in_class,
                          n_hits,
                          lower.tail=FALSE)
            
            # Add metrics to the result data frame.
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
    
    # Adjust p-values.
    results$AdjPVal = p.adjust(results$PVal, method="bonferroni")
    
    return(results)
}

# Write output table.
dir.create("output", showWarnings=FALSE)

# Perform enrichment only for terms present more than 10 times.
# Terms which are barely representated will give anecdotal results at best,
# and their large number will cause a large increase in p-values at the 
# multiple-tests correction step.
large_drug_classes = names(table(all_data$Class)[table(all_data$Class) >= 10])

# Start by testing everything at once and applying the appropriate 
# FDR correction.
results = test_drug_classes(all_data, large_drug_classes, hit_classes)
write.table(results[order(results$PVal),], file="output/Drug class enrichment (large classes).txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

# Perform enrichments one class at a time to see how corrected p-values are affected.
for(hit_class in names(hit_classes)) {
    results = test_drug_classes(all_data, large_drug_classes, hit_classes[hit_class])
    write.table(results[order(results$PVal),], sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE, 
                file=paste0("output/Drug class enrichment - ", hit_class, " - (large classes).txt"))
}


# Perform enrichment of terms present less than 10 times
# to identify anecdotal hits in case they may be cross-referenced
# with information from the litterature later on.
small_drug_classes = names(table(all_data$Class)[table(all_data$Class) < 10])
results = test_drug_classes(all_data, small_drug_classes, hit_classes)
write.table(results[order(results$PVal),], file="output/Drug class enrichment (small classes).txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
