library(phyloseq)
library(mirlyn)
data(example)

test_that("mirl", {
  mirlobj <- mirl(example, rep = 5, set.seed = 1)
  expect_type(mirlobj, "list")
})

test_that("alphacone", {
  expect_type(alphacone(example, rep = 5, steps = c(0.25, 0.5, 0.75)), "list")
})

test_that("alphawhich", {
  adf <- alphadivDF(mirlobj)
  expect_type(alphawichVis(adf, "time"), "list")
})

test_that("bartax", {
  cols <- c("black","darkgoldenrod1","dodgerblue","deeppink4","chartreuse3",
     "burlywood4","navy","blueviolet", "tan2","lavenderblush3","cyan4")
  expect_type(bartax(example, "Sample", lowabun = 0.01, yvar = "Abundance", taxrank = "Phylum", cols), "list")
})

test_that("betamat", {
  pcaobj <- betamatPCA(mirlobj, dsim = "jaccard")
  expect_type(pcaobj, "list")
  expect_type(betamatPCAvis(pcaobj, groups = LETTERS[1:6], colours = c("#000000", "#E69F00", "#0072B2", "#009E73", "#F0E442", "#D55E00")), "list")
})

# test_that("divana", {
#   expect_type(divana(example, reps = 3), "list")
# })

test_that("rarecurve", {
  rwre <- rarefy_whole_rep(example, rep = 5, steps = c(0.25, 0.5, 0.75))
  expect_type(rarecurve(rwre), "list")
})
