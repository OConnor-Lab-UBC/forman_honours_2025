
weights <- read.csv("/Users/ericaforman/Desktop/anemone_aggregation/processed data - rawdata (3).csv")

counts <- read.csv("/Users/ericaforman/Desktop/anemone_aggregation/counts.csv")

field <- read.csv("/Users/ericaforman/Desktop/anemone_aggregation/raw collections data - field data.csv")

field$anemone_weight <- as.numeric(field$anemone_weight)
class(field$anemone_weight)


#plot anemone wet weight vs aggregate size from field data
field_grouped <- field|>
  group_by(quadrat, aggregate_size, site)|>
  summarize(mean_weight = mean(anemone_weight))

ggplot(field_grouped, aes(x=aggregate_size, y=mean_weight, colour = site))+
  geom_point()

anemone_wet_weight_aggregate_model <- lmer(mean_weight ~ aggregate_size +
                                             (1 | site), data = field_grouped)
summary(anemone_wet_weight_aggregate_model)

install.packages("simr")
library(simr)

extended_model <- simr::extend(anemone_wet_weight_aggregate_model, 
                               along = "site",
                               n = 1000)

power_curve <- powerCurve(extended_model,
                          test=simr::fixed("aggregate_size", method="lr"),
                          along = "site",
                          breaks = c(0, 200, 400, 600, 800, 1000),
                          nsim = 100,
                          alpha=.05,
                          progress=TRUE)

summary(power_curve)


#plot algae weight vs anemone weight 
grouped_weights <- weights|>
group_by(sample_number, sample_identity, site)|>
  summarise(total_empty = sum(emptyboat_weight), 
            total_48h = sum(X48h_weight), 
            total_72h = sum(X72h_weight))|>
  mutate(sample_weight_72h=(total_72h-total_empty))

wide_weights <- grouped_weights|>
  pivot_wider(names_from = sample_identity, 
              values_from = sample_weight_72h,
              id_cols = c(site, sample_number))

ggplot(wide_weights, aes(x=anemone, y=algae, colour=site))+
  geom_point()

algae_weight_anemone_weight_model <- lmer(algae ~ anemone + (1|site), data = wide_weights)
summary(algae_weight_anemone_weight_model)

extended_model <- simr::extend(algae_weight_anemone_weight_model, 
                               along = "site",
                               n = 200)

power_curve <- powerCurve(extended_model,
         test=simr::fixed("anemone", method="lr"),
         along = "site",
         breaks = c(0, 50, 100, 150, 200),
         nsim = 100,
         alpha=.05,
         progress=TRUE)

summary(power_curve)


#clean counts data to get actual counts 
mutated_counts <- counts|>
  mutate(dilution_factor=subsample_total_v/subsample_algae_v)|>
  group_by(sample_number, dilution_factor, subsample_algae_v, subsample_total_v)|>
  summarise(avg_count=mean(average))|>
  replace_na(list(subsample_total_v = 0, subsample_algae_v = 0, dilution_factor = 1))|>
  mutate(total_count=dilution_factor*avg_count*(3000-subsample_algae_v))

#plot counts vs algae weight 
weights_and_counts <- mutated_counts|>
  full_join(wide_weights, by = "sample_number")

ggplot(weights_and_counts, aes(x=algae, y=total_count, colour=site))+
  geom_point()

algae_counts_algae_weight_model <- lmer(total_count ~ algae + (1|site), 
                                        data = weights_and_counts)
summary(algae_counts_algae_weight_model)

#plot counts vs anemone weight
ggplot(weights_and_counts, aes(x=anemone, y=total_count, colour=site))+
  geom_point()

algae_counts_anemone_weight_model <- lmer(total_count ~ anemone + (1|site),
                                          data = weights_and_counts)
summary(algae_counts_anemone_weight_model)

#plot algae weight vs aggregate size 
field_and_weights <- wide_weights|>
  select(-site) |>
  left_join(field, by = c("sample_number"))

ggplot(field_and_weights, aes(x=aggregate_size, y=algae, colour = site))+
  geom_point()

algae_weight_aggregate_size_model <- lmer(algae ~ aggregate_size + (1|site),
                                          data = field_and_weights)

#plot anemone weight vs aggregate size
ggplot(field_and_weights, aes(x=aggregate_size, y=anemone, colour = site))+
  geom_point()

  