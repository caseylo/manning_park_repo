library(ggplot2)


model_data_late <- read_csv("manning_park_data/model_data.csv")

cor.test(model_data_not07$Year, model_data_not07$Elevation)

model_data_not07 <- model_data %>%
  filter(Year < 2025) %>%
  filter(Year != 2007)

ggplot(model_data, aes(x = Year, y = Elevation, colour = ProjectID)) +
   geom_point()

model_data_late <- model_data %>%
  filter(Year > 1993)

mod8_full <- glmer(
  present ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + Lat_sc + Lon_sc,
  data = model_data_late,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)


summary(mod8_full)
ranef(mod8_full)