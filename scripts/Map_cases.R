library(dplyr)
library(tibble)
library(ggplot2)
library(sf)
library(maptiles)
library(tidyterra)

# Static map of stillbirths and neonatal/infant deaths against all births.
# Needs short_raw, births_wash and analysis_df (from scripts/Load_data.R).
# Produces: case_map (ggplot), map_summary (placement counts), map_points (plotted coordinates)
#
# Placement rules (GPS is missing for most records):
#   1. Use the recorded GPS (start-of-survey, else consent, else end geopoint) if it is inside
#      the study area, within max_km of its subcounty centre and inside its subcounty polygon.
#   2. Otherwise place at the median GPS of the village (else parish, else subcounty),
#      with +/- jitter degrees of uniform noise (about 300 m). Centroids use valid GPS only.
# Subcounty polygons: geoBoundaries ADM4 (2019). The recorded subcounties follow the newer
# structure, so the Jinja divisions are unions of old units chosen from where the GPS points fall
# (approximate), and Namagera Town Council has no polygon (its GPS is not polygon-tested).

max_km <- 20
jit    <- 0.003
seed   <- 4817
geojson_url <- "https://github.com/wmgeolab/geoBoundaries/raw/9469f09/releaseData/gbOpen/UGA/ADM4/geoBoundaries-UGA-ADM4_simplified.geojson"
geojson_file <- file.path("output", "geoBoundaries-UGA-ADM4_simplified.geojson")

old_s2 <- sf_use_s2(FALSE)

ok  <- \(lat, lon) !is.na(lat) & lat > -0.5 & lat < 1.5 & lon > 33 & lon < 34.2
hav <- \(lat1, lon1, lat2, lon2) {
  r <- pi / 180
  a <- sin((lat2 - lat1) * r / 2)^2 + cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  12742 * asin(sqrt(a))
}

# One row per participant: admin units and best available GPS
loc <- short_raw |>
  transmute(
    survey_participant_number,
    district = survey_participant_district, subcounty = survey_participant_subcounty, parish, village,
    lat = coalesce(survey_start_geopoint_latitude, consent_form_consent_geopoint_latitude, end_geopoint_latitude),
    lon = coalesce(survey_start_geopoint_longitude, consent_form_consent_geopoint_longitude, end_geopoint_longitude)
  ) |>
  distinct(survey_participant_number, .keep_all = TRUE)

# Subcounty polygons ---------------------------------------------------------
if (!file.exists(geojson_file)) {
  dir.create(dirname(geojson_file), showWarnings = FALSE, recursive = TRUE)
  download.file(geojson_url, geojson_file, quiet = TRUE, mode = "wb")
}
adm4 <- st_read(geojson_file, quiet = TRUE) |> st_make_valid() |> mutate(poly_id = row_number())

xwalk <- tribble(
  ~subcounty, ~shapeName, ~poly_id,
  "BUSEDDE", "Busedde", 1201, "BUTAGAYA", "Butagaya", 1170, "BUWENGE", "Buwenge", 1171,
  "BUWENGE TOWN COUNCIL", "Buwenge Town Council", 1172, "BUYENGO TOWN COUNCIL", "Buyengo", 1199,
  "KAKIRA TOWN COUNCIL", "Kakira Town Council", 1167,
  "JINJA NORTH DIVISION", "Mafubira", 1168, "JINJA NORTH DIVISION", "Budondo", 1169,
  "JINJA NORTH DIVISION", "Bugembe Town Council", 1166,
  "JINJA SOUTH DIVISION", "Mpumudde Division", 1174, "JINJA SOUTH DIVISION", "Central Division", 1162,
  "JINJA SOUTH DIVISION", "Masese Walukuba", 1173
)
stopifnot(all(adm4$shapeName[xwalk$poly_id] == xwalk$shapeName))

bounds <- adm4 |>
  inner_join(xwalk |> select(subcounty, poly_id), by = "poly_id") |>
  group_by(subcounty) |>
  summarise(geometry = st_union(geometry), .groups = "drop") |>
  suppressWarnings()

# Valid GPS: in study area, near its subcounty centre, inside its polygon ---
gps <- loc |> filter(ok(lat, lon))
sub_centre <- gps |> group_by(district, subcounty) |>
  summarise(slat = median(lat), slon = median(lon), .groups = "drop")
gps <- gps |>
  left_join(sub_centre, by = c("district", "subcounty")) |>
  filter(hav(lat, lon, slat, slon) <= max_km)

in_poly <- gps |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |>
  select(survey_participant_number, subcounty) |>
  group_split(subcounty) |>
  purrr::map(\(d) {
    poly <- bounds[bounds$subcounty == d$subcounty[1], ]
    tibble(survey_participant_number = d$survey_participant_number,
           outside = if (nrow(poly)) lengths(suppressMessages(st_intersects(d, poly))) == 0 else FALSE)
  }) |>
  bind_rows()
gps <- gps |>
  left_join(in_poly, by = "survey_participant_number") |>
  filter(!outside)

# Centroids from valid GPS only ----------------------------------------------
cent <- \(by) gps |> group_by(across(all_of(by))) |>
  summarise(clat = median(lat), clon = median(lon), .groups = "drop")
c_vil <- cent(c("district", "subcounty", "parish", "village")) |> rename(vlat = clat, vlon = clon)
c_par <- cent(c("district", "subcounty", "parish")) |> rename(plat = clat, plon = clon)
c_sub <- cent(c("district", "subcounty")) |> rename(slat = clat, slon = clon)
valid_ids <- gps$survey_participant_number

place <- function(ids) {
  set.seed(seed)
  tibble(survey_participant_number = ids) |>
    left_join(loc, by = "survey_participant_number") |>
    left_join(c_vil, by = c("district", "subcounty", "parish", "village")) |>
    left_join(c_par, by = c("district", "subcounty", "parish")) |>
    left_join(c_sub, by = c("district", "subcounty")) |>
    mutate(
      level = case_when(
        survey_participant_number %in% valid_ids ~ "GPS",
        !is.na(vlat) ~ "village", !is.na(plat) ~ "parish", !is.na(slat) ~ "subcounty", TRUE ~ "none"),
      lat2 = case_when(level == "GPS" ~ lat, level == "village" ~ vlat, level == "parish" ~ plat, level == "subcounty" ~ slat),
      lon2 = case_when(level == "GPS" ~ lon, level == "village" ~ vlon, level == "parish" ~ plon, level == "subcounty" ~ slon),
      lat2 = lat2 + if_else(level == "GPS", 0, runif(n(), -jit, jit)),
      lon2 = lon2 + if_else(level == "GPS", 0, runif(n(), -jit, jit))
    ) |>
    filter(level != "none") |>
    select(survey_participant_number, level, lat2, lon2)
}

# Points to draw ---------------------------------------------------------------
cases <- analysis_df |>
  filter(stillbirth == 1 | died_any) |>
  transmute(survey_participant_number, outcome = if_else(stillbirth == 1, "Stillbirth", "Death")) |>
  inner_join(place(analysis_df$survey_participant_number[analysis_df$stillbirth == 1 | analysis_df$died_any]),
             by = "survey_participant_number")
others <- \(ids) place(setdiff(ids, cases$survey_participant_number))

map_points <- bind_rows(
  mutate(others(births_wash$survey_participant_number), group = "Births"),
  mutate(others(analysis_df$survey_participant_number), group = "Births with 28-day follow-up"),
  mutate(cases, group = outcome)
)
map_summary <- count(cases, level, outcome)

# Plot -------------------------------------------------------------------------
to_sf <- \(d) st_as_sf(d, coords = c("lon2", "lat2"), crs = 4326) |> st_transform(3857)
bounds_m <- st_transform(bounds, 3857)
bb <- st_bbox(bounds_m)
tiles <- get_tiles(st_buffer(st_as_sfc(st_bbox(bounds)), 0.02) |> suppressWarnings(),
                   provider = "Esri.WorldGrayCanvas", zoom = 11, crop = TRUE, project = FALSE)

case_map <- ggplot() +
  geom_spatraster_rgb(data = tiles, maxcell = 5e6) +
  geom_sf(data = bounds_m, fill = NA, color = "grey25", linewidth = 0.3) +
  geom_sf(data = to_sf(filter(map_points, group == "Births")),
          color = "grey40", alpha = 0.15, size = 0.6, shape = 16) +
  geom_sf(data = to_sf(filter(map_points, group == "Births with 28-day follow-up")),
          color = "steelblue", alpha = 0.2, size = 0.6, shape = 16) +
  geom_sf(data = to_sf(cases), aes(color = outcome), size = 1.8, alpha = 0.9) +
  scale_color_manual(values = c(Death = "firebrick", Stillbirth = "navy"), name = NULL) +
  coord_sf(xlim = bb[c("xmin", "xmax")], ylim = bb[c("ymin", "ymax")], expand = TRUE) +
  theme_void() +
  theme(legend.position = "inside", legend.position.inside = c(0.12, 0.12))

sf_use_s2(old_s2)
