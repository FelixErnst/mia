test_that("getMediation", {

  skip_if_not(require("miaTime", quietly = TRUE))
  data("hitchip1006", package = "miaTime")
  tse <- hitchip1006

  tse <- agglomerateByRank(tse, rank = "Phylum")
  tse$bmi_group <- as.numeric(tse$bmi_group)

  ### Batch 1: check errors when missing or wrong arguments ###
  expect_error(
    getMediation(tse, outcome = "bmi_group", treatment = "nationality",
                 mediator = "diversity", assay.type = "counts"),
    "The arguments mediator, assay.type and dimred are mutually exclusive, but 2 were provided."
  )

  expect_error(
    getMediation(tse, outcome = "wrong_name", treatment = "nationality", mediator = "diversity"),
    "wrong_name not found in colData(x).", fixed = TRUE
  )

  expect_error(
    getMediation(tse, outcome = "bmi_group", treatment = "wrong_name", mediator = "diversity"),
    "wrong_name not found in colData(x).", fixed = TRUE
  )

  expect_error(
    getMediation(tse, outcome = "bmi_group", treatment = "nationality", mediator = "wrong_name"),
    "wrong_name not found in colData(x).", fixed = TRUE
  )

  expect_error(
    getMediation(tse, outcome = "bmi_group", treatment = "nationality", assay.type = "wrong_name"),
    "'assay.type' must be a valid name of assays(x)", fixed = TRUE
  )

  expect_error(
    getMediation(tse, outcome = "bmi_group", treatment = "nationality", dimred = "wrong_name"),
    "wrong_name not found in reducedDims(x).", fixed = TRUE
  )

  expect_error(
    getMediation(tse, outcome = "bmi_group", treatment = "nationality",
                 mediator = "diversity", boot = TRUE, sims = 1),
    "Too many treatment levels. Consider specifing a treat.value and a control.value"
  )

  expect_error(
    getMediation(tse, outcome = "bmi_group", treatment = "nationality", mediator = "diversity",
                 boot = TRUE, sims = 1, treat.value = "wrong_value", control.value = "wrong_value"),
    "treat.value and/or control.value not found in the levels of the treatment variable."
  )

  skip_if_not_installed("mediation")
  ### Batch 2: check equality between base function and wrapper ###
  set.seed(123)
  med_df <- getMediation(tse, outcome = "bmi_group", treatment = "nationality", mediator = "diversity",
                         treat.value = "Scandinavia", control.value = "CentralEurope",
                         boot = TRUE, sims = 1, add.metadata = TRUE)

  df <- data.frame(Outcome = tse$bmi_group, Treatment = tse$nationality, Mediator = tse$diversity)
  df <- na.omit(df)
  df <- df[df$Treatment %in% c("Scandinavia", "CentralEurope"), ]

  fit_m <- lm(Mediator ~ Treatment, data = df)
  fit_dv <- glm(Outcome ~ Treatment + Mediator, data = df)

  set.seed(123)
  med_out <- mediate(fit_m, fit_dv, treat = "Treatment", mediator = "Mediator",
                     treat.value = "Scandinavia", control.value = "CentralEurope",
                     boot = TRUE, sims = 1)

  expect_equal(attr(med_df, "model_metadata")[[1]][["d.avg"]], med_out$d.avg)
  expect_equal(attr(med_df, "model_metadata")[[1]][["d.avg.p"]], med_out$d.avg.p)
  expect_equal(attr(med_df, "model_metadata")[[1]][["z.avg"]], med_out$z.avg)
  expect_equal(attr(med_df, "model_metadata")[[1]][["z.avg.p"]], med_out$z.avg.p)

  ### Batch 3: check output format and dimensionality with respect to SE ###
  med_df <- getMediation(tse, outcome = "bmi_group", treatment = "nationality", assay.type = "counts",
                         treat.value = "Scandinavia", control.value = "CentralEurope",
                         boot = TRUE, sims = 1)

  expect_named(tse, med_df[["mediator"]])

  expect_named(med_df, c("mediator", "acme", "acme_pval", "acme_lower",
      "acme_upper", "ade", "ade_pval", "ade_lower", "ade_upper", "total",
      "total_lower", "total_upper", "total_pval", "acme_padj", "ade_padj",
      "total_padj"))
})
