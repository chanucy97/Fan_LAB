library(png)
library(grid)
library(ragg)
library(svglite)

root <- normalizePath(Sys.getenv("FANLAB_RESEARCH2_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "final_manuscript_planning_20260605", "figures_redrawn", "fig1_cohort_design_overview")

panel_a_path <- file.path(out_dir, "Fig1_cohort_design_overview.png")
panel_b_path <- file.path(out_dir, "Fig1B_coverage_heatmap.png")

panel_a <- readPNG(panel_a_path)
panel_b <- readPNG(panel_b_path)

width_mm <- 183
gap_mm <- 4
a_height_mm <- width_mm * dim(panel_a)[1] / dim(panel_a)[2]
b_height_mm <- width_mm * dim(panel_b)[1] / dim(panel_b)[2]
height_mm <- a_height_mm + gap_mm + b_height_mm

draw_fig1 <- function() {
  grid.newpage()
  pushViewport(viewport(width = unit(1, "npc"), height = unit(1, "npc")))
  grid.rect(gp = gpar(fill = "white", col = NA))

  a_h <- a_height_mm / height_mm
  b_h <- b_height_mm / height_mm
  gap_h <- gap_mm / height_mm

  pushViewport(viewport(
    x = unit(0.5, "npc"),
    y = unit(1 - a_h / 2, "npc"),
    width = unit(1, "npc"),
    height = unit(a_h, "npc"),
    clip = "on"
  ))
  grid.raster(panel_a, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)
  grid.text(
    "A",
    x = unit(5, "mm"),
    y = unit(1, "npc") - unit(5, "mm"),
    just = c("left", "top"),
    gp = gpar(fontfamily = "Arial", fontsize = 16, fontface = "bold", col = "black")
  )
  popViewport()

  pushViewport(viewport(
    x = unit(0.5, "npc"),
    y = unit(b_h / 2, "npc"),
    width = unit(1, "npc"),
    height = unit(b_h, "npc"),
    clip = "on"
  ))
  grid.raster(panel_b, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)
  popViewport()

  popViewport()
}

base <- file.path(out_dir, "Fig1_complete_cohort_coverage")

ragg::agg_png(paste0(base, ".png"), width = width_mm / 25.4, height = height_mm / 25.4, units = "in", res = 600, background = "white")
draw_fig1()
dev.off()

ragg::agg_tiff(paste0(base, ".tiff"), width = width_mm / 25.4, height = height_mm / 25.4, units = "in", res = 600, background = "white")
draw_fig1()
dev.off()

grDevices::cairo_pdf(paste0(base, ".pdf"), width = width_mm / 25.4, height = height_mm / 25.4, family = "Arial", bg = "white")
draw_fig1()
dev.off()

svglite::svglite(paste0(base, ".svg"), width = width_mm / 25.4, height = height_mm / 25.4, bg = "white")
draw_fig1()
dev.off()

writeLines(
  c(
    "Figure 1 complete assembly",
    paste0("Panel A: ", panel_a_path),
    paste0("Panel B: ", panel_b_path),
    paste0("Width mm: ", round(width_mm, 2)),
    paste0("Height mm: ", round(height_mm, 2)),
    paste0("Panel A height mm: ", round(a_height_mm, 2)),
    paste0("Panel B height mm: ", round(b_height_mm, 2))
  ),
  con = file.path(out_dir, "Fig1_complete_assembly_notes.txt")
)

message("Wrote complete Figure 1 to: ", out_dir)
