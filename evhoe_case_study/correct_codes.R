## Fit a spatial model on Hake with sdmTMB
##----------------------------------------
# To run the codes, you'll need to have run the codes in cours.qmd
# (at least the chunk for loading the data and build `Catch_sf_species`)

## Load packages
library(sdmTMB)
library(sf)

## Build data frame to fit sdmTMB models
train_df <- cbind(Catch_sf_species,st_coordinates(Catch_sf_species))
train_df_2 <- train_df |> 
  filter(!is.na(Depth)) |> # quick fix : filter few NA values in Depth
  as.data.frame() |> 
  dplyr::select(TotalBiom,Year,Depth,X,Y)

## Make mesh
mesh <- make_mesh(train_df_2, xy_cols = c("X", "Y"), n_knots = 100)
plot(mesh)

## Fit spatial model
fit <- sdmTMB(
  TotalBiom ~ s(Depth),
  data = train_df_2,
  mesh = mesh,
  family = tweedie(link = "log"),
  spatial = "on"
)

## Check that the model has converged
sanity(fit)

## Look at the depth effect
ggeffects::ggpredict(fit, "Depth") |> plot()

## Create interpolation grid
x_range <- range(train_df$X)
y_range <- range(train_df$Y)
res <- 0.05
grid <- expand.grid(
  X = seq(from = x_range[1], to = x_range[2], by = res),
  Y = seq(from = y_range[1], to = y_range[2], by = res)
)
point_grid_sf <- st_as_sf(
  grid,
  coords = c("X", "Y"),
  crs = st_crs(mapBase)
)

# Now filter the locations that are within the area of interest
pts_sf <- st_as_sf(
  train_df_2,
  coords = c("X", "Y"),
  crs = st_crs(mapBase)
)
hull <- concaveman( # build buffer around data points
  pts_sf,
  concavity = 2,      # higher = smoother
  length_threshold = 0
)

point_grid_sf_2 <- st_intersection(point_grid_sf,hull)
grid_2 <- st_coordinates(point_grid_sf_2) |> data.frame()

## Create new data for prediction
newdata <- data.frame(Depth = rep(100,nrow(grid_2)), # we fix the depth at this stage
  X = grid_2$X,
  Y = grid_2$Y)

## Predict and plot
p <- predict(fit, newdata = newdata)
p_sf <- st_as_sf(p,coords = c("X","Y"),crs = st_crs(mapBase))
plot_without_depth <- ggplot() +
  geom_sf(data = p_sf, aes(col = exp(est)))+
  scale_color_distiller(palette = "Spectral")+
  geom_sf(data = mapBase) +
  # geom_sf(data = Catch_sf_species,aes(size = TotalBiom),col="grey",alpha=0.05)+
  coord_sf(xlim = xlims, ylim = ylims, expand = FALSE)+
  ggtitle("Predictions without the depth effect")

## But big problem --> no interpolated values of depth
## So we predict depth on the interpolated grid for a model for depth
fit_model_depth <- sdmTMB(
  Depth ~ 1,
  data = train_df_2[,c("Depth","X","Y")],
  mesh = mesh,
  family = Gamma(link = "inverse"),
  spatial = "on",
)
sanity(fit_model_depth)

## Predict and plot depth
newdata_pred_depth <- data.frame(X = grid_2$X,
  Y = grid_2$Y)
p_depth <- predict(fit_model_depth, newdata = newdata_pred_depth)
p_depth_sf <- st_as_sf(p_depth,coords = c("X","Y"),crs = st_crs(mapBase))
plot_depth <- ggplot() +
  geom_sf(data = p_depth_sf, aes(col = 1/est))+
  scale_color_distiller(palette = "Spectral",limits = c(0,550))+
  geom_sf(data = mapBase) +
  # geom_sf(data = Catch_sf_species,aes(size = Depth),col="grey",alpha=0.05)+
  coord_sf(xlim = xlims, ylim = ylims, expand = FALSE)+
  ggtitle("Predicted values od depth")

## Now we add depth to the predict dataset
newdata_with_depth <- data.frame(Depth = 1/p_depth$est,
  X = grid_2$X,
  Y = grid_2$Y)

## Predict and plot
p_with_depth <- predict(fit, newdata = newdata_with_depth)
p_with_depth_sf <- st_as_sf(p_with_depth,coords = c("X","Y"),crs = st_crs(mapBase))
plot_with_depth <- ggplot() +
  geom_sf(data = p_with_depth_sf, aes(col = exp(est)))+
  scale_color_distiller(palette = "Spectral")+
  geom_sf(data = mapBase) +
  # geom_sf(data = Catch_sf_species,aes(size = TotalBiom),col="grey",alpha=0.25)+
  coord_sf(xlim = xlims, ylim = ylims, expand = FALSE)+
  ggtitle("Predictions with the depth effect")

cowplot::plot_grid(plot_without_depth,plot_with_depth,plot_depth,nrow = 1)

## TO DO :)
## Add a spatio-temporal term (the data are spatio-temporal, but at this point the model is only spatial)
## Check the residuals --> look at https://sdmtmb.github.io/sdmTMB/articles/residual-checking.html
## Asssess sensitivity to the mesh
## Apply it to other species