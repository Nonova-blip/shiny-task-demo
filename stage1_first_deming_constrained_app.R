# Stage 1-first Deming-constrained signal task prototype
# ------------------------------------------------------------
# Purpose:
#   A simpler two-stage Shiny prototype for the dot-cloud task.
#
# Design logic:
#   1. Generate Stage 1 first around a latent target line.
#   2. Accept/reject Stage 1 using its Deming/orthogonal slope.
#   3. Add independent Stage 2 dots from the same latent line.
#   4. Accept/reject the full Stage 2 display using its Deming/orthogonal slope.
#
# Why this version exists:
#   - Avoids the earlier Stage 2-first subset-selection artifact.
#   - Restores Deming-based regeneration so realized displayed slopes do not drift too far.
#   - Keeps the app minimal: no OLS lines, no extra participant-view machinery.
#
# To run:
#   install.packages("shiny")  # if needed
#   shiny::runApp("stage1_first_deming_constrained_app.R")

library(shiny)

# ---------- Deming / orthogonal regression ----------
deming_fit <- function(x, y, lambda = 1) {
  xbar <- mean(x)
  ybar <- mean(y)

  Sxx <- mean((x - xbar)^2)
  Syy <- mean((y - ybar)^2)
  Sxy <- mean((x - xbar) * (y - ybar))

  if (abs(Sxy) < 1e-12) {
    slope <- 0
  } else {
    slope <- (Syy - lambda * Sxx +
                sqrt((Syy - lambda * Sxx)^2 + 4 * lambda * Sxy^2)) /
      (2 * Sxy)
  }

  intercept <- ybar - slope * xbar
  list(intercept = intercept, slope = slope)
}

# ---------- Draw bounded orthogonal noise ----------
draw_bounded_noise <- function(n, noise_sd, bound, max_iter = 5000) {
  for (i in seq_len(max_iter)) {
    d <- rnorm(n, mean = 0, sd = noise_sd)
    if (all(abs(d) <= bound)) return(d)
  }
  NULL
}

# ---------- Generate dots around a latent line using orthogonal noise ----------
generate_cloud <- function(target_slope,
                           n_dots,
                           noise_sd,
                           bound,
                           t_range = 0.42,
                           center_x = 0.5,
                           center_y = 0.5,
                           source = "dots") {

  # Even spacing along the latent line keeps the intended structure more stable.
  # Randomness enters through orthogonal noise, not through uneven x coverage.
  t_vals <- seq(-t_range, t_range, length.out = n_dots)

  v_line <- c(1, target_slope)
  v_line <- v_line / sqrt(sum(v_line^2))

  v_perp <- c(-target_slope, 1)
  v_perp <- v_perp / sqrt(sum(v_perp^2))

  d <- draw_bounded_noise(n_dots, noise_sd, bound)
  if (is.null(d)) return(NULL)

  x <- center_x + t_vals * v_line[1] + d * v_perp[1]
  y <- center_y + t_vals * v_line[2] + d * v_perp[2]

  # Keep displays inside the plotting frame. This is a pragmatic prototype constraint.
  if (any(x < 0 | x > 1 | y < 0 | y > 1)) return(NULL)

  out <- data.frame(
    x = x,
    y = y,
    t = t_vals,
    orthogonal_noise = d,
    source = source,
    stringsAsFactors = FALSE
  )

  dem <- deming_fit(out$x, out$y)
  attr(out, "deming_slope") <- dem$slope
  attr(out, "deming_intercept") <- dem$intercept
  out
}

# ---------- Main two-stage generator ----------
make_two_stage_trial <- function(target_slope = 0.10,
                                 n_stage1 = 8,
                                 n_stage2 = 20,
                                 noise_sd_stage1 = 0.12,
                                 noise_sd_added = 0.12,
                                 bound = 0.28,
                                 stage1_tol = 0.08,
                                 stage2_tol = 0.03,
                                 require_stage2_not_worse = TRUE,
                                 max_stage1_tries = 2000,
                                 max_added_tries = 2000) {

  if (n_stage1 >= n_stage2) stop("Stage 1 dots must be fewer than Stage 2 dots.")
  n_added <- n_stage2 - n_stage1

  for (s1_try in seq_len(max_stage1_tries)) {
    stage1 <- generate_cloud(
      target_slope = target_slope,
      n_dots = n_stage1,
      noise_sd = noise_sd_stage1,
      bound = bound,
      source = "Stage 1"
    )

    if (is.null(stage1)) next

    dem1_slope <- attr(stage1, "deming_slope")
    s1_error <- abs(dem1_slope - target_slope)

    # Stage 1 is constrained, but usually more loosely than Stage 2.
    # This avoids wildly off-slope Stage 1 trials without forcing Stage 1 to be perfect.
    if (s1_error > stage1_tol) next

    for (add_try in seq_len(max_added_tries)) {
      added <- generate_cloud(
        target_slope = target_slope,
        n_dots = n_added,
        noise_sd = noise_sd_added,
        bound = bound,
        source = "Added at Stage 2"
      )

      if (is.null(added)) next

      stage2 <- rbind(stage1, added)
      dem2 <- deming_fit(stage2$x, stage2$y)
      dem2_slope <- dem2$slope
      s2_error <- abs(dem2_slope - target_slope)

      if (s2_error > stage2_tol) next
      if (require_stage2_not_worse && s2_error > s1_error) next

      attr(stage1, "target_slope") <- target_slope
      attr(stage1, "s1_try") <- s1_try
      attr(stage1, "stage1_error") <- s1_error

      attr(stage2, "target_slope") <- target_slope
      attr(stage2, "deming_slope") <- dem2_slope
      attr(stage2, "deming_intercept") <- dem2$intercept
      attr(stage2, "stage2_error") <- s2_error
      attr(stage2, "s1_try") <- s1_try
      attr(stage2, "added_try") <- add_try
      attr(stage2, "n_added") <- n_added

      return(list(stage1 = stage1, stage2 = stage2, added = added))
    }
  }

  NULL
}

# ---------- Plot helper ----------
draw_trial_panel <- function(d,
                             title_text,
                             target_slope,
                             show_target = TRUE,
                             show_deming = TRUE,
                             highlight_added = FALSE) {
  point_cols <- rep("black", nrow(d))
  if (highlight_added && "source" %in% names(d)) {
    point_cols[d$source == "Added at Stage 2"] <- "gray70"
  }

  plot(
    d$x, d$y,
    xlim = c(0, 1), ylim = c(0, 1),
    pch = 19, cex = 1.65,
    col = point_cols,
    xlab = "", ylab = "",
    axes = FALSE,
    main = title_text
  )
  box()

  if (show_target) {
    target_intercept <- 0.5 - target_slope * 0.5
    abline(a = target_intercept, b = target_slope, col = "blue", lty = 2, lwd = 2.5)
  }

  if (show_deming) {
    abline(a = attr(d, "deming_intercept"),
           b = attr(d, "deming_slope"),
           col = "darkgreen", lwd = 2.5)
  }
}

# ---------- UI ----------
ui <- fluidPage(
  titlePanel("Two-Stage Dot Task: Stage 1 First + Deming Constraints"),

  sidebarLayout(
    sidebarPanel(
      tags$h4("Latent structure"),
      sliderInput("slope", "Target slope:",
                  min = -0.40, max = 0.40, value = 0.10, step = 0.01),

      tags$h4("Dots"),
      sliderInput("n1", "Stage 1 dots:",
                  min = 4, max = 25, value = 8, step = 1),
      sliderInput("n2", "Stage 2 total dots:",
                  min = 6, max = 60, value = 20, step = 1),

      tags$h4("Noise"),
      sliderInput("noise1", "Stage 1 orthogonal noise SD:",
                  min = 0.01, max = 0.30, value = 0.12, step = 0.01),
      sliderInput("noise_added", "Added-dot orthogonal noise SD:",
                  min = 0.01, max = 0.30, value = 0.12, step = 0.01),
      sliderInput("bound", "Orthogonal noise bound:",
                  min = 0.03, max = 0.45, value = 0.28, step = 0.01),

      tags$h4("Deming regeneration rules"),
      sliderInput("tol1", "Stage 1 slope tolerance:",
                  min = 0.005, max = 0.20, value = 0.08, step = 0.005),
      sliderInput("tol2", "Stage 2 slope tolerance:",
                  min = 0.005, max = 0.10, value = 0.03, step = 0.005),
      checkboxInput("not_worse", "Require Stage 2 slope to be at least as close to target as Stage 1", TRUE),

      tags$h4("Display"),
      checkboxInput("show_target", "Show latent target line", TRUE),
      checkboxInput("show_deming", "Show Deming fit", TRUE),
      checkboxInput("highlight_added", "Show added Stage 2 dots in gray", TRUE),

      numericInput("tries1", "Max Stage 1 tries:", value = 2000, min = 100, step = 100),
      numericInput("tries2", "Max added-dot tries:", value = 2000, min = 100, step = 100),

      actionButton("regen", "Generate Trial")
    ),

    mainPanel(
      plotOutput("stage1plot", height = "330px"),
      plotOutput("stage2plot", height = "330px"),
      tags$hr(),
      verbatimTextOutput("info")
    )
  )
)

# ---------- Server ----------
server <- function(input, output, session) {

  trial <- eventReactive(input$regen, {
    validate(need(input$n1 < input$n2, "Stage 1 dots must be fewer than Stage 2 dots."))

    make_two_stage_trial(
      target_slope = input$slope,
      n_stage1 = input$n1,
      n_stage2 = input$n2,
      noise_sd_stage1 = input$noise1,
      noise_sd_added = input$noise_added,
      bound = input$bound,
      stage1_tol = input$tol1,
      stage2_tol = input$tol2,
      require_stage2_not_worse = input$not_worse,
      max_stage1_tries = input$tries1,
      max_added_tries = input$tries2
    )
  }, ignoreInit = FALSE)

  output$stage1plot <- renderPlot({
    tr <- trial()
    validate(need(!is.null(tr), "No valid trial found. Loosen tolerances, reduce noise, increase bound, or increase tries."))

    draw_trial_panel(
      d = tr$stage1,
      title_text = "Stage 1: initial evidence",
      target_slope = input$slope,
      show_target = input$show_target,
      show_deming = input$show_deming,
      highlight_added = FALSE
    )
  })

  output$stage2plot <- renderPlot({
    tr <- trial()
    validate(need(!is.null(tr), "No valid trial found. Loosen tolerances, reduce noise, increase bound, or increase tries."))

    draw_trial_panel(
      d = tr$stage2,
      title_text = "Stage 2: Stage 1 plus added evidence",
      target_slope = input$slope,
      show_target = input$show_target,
      show_deming = input$show_deming,
      highlight_added = input$highlight_added
    )
  })

  output$info <- renderPrint({
    tr <- trial()

    if (is.null(tr)) {
      cat("No valid trial found.\n")
      cat("Try: loosen Stage 1/Stage 2 tolerances, lower noise, increase bound, or increase tries.\n")
      return()
    }

    d1 <- tr$stage1
    d2 <- tr$stage2

    cat("Generation logic: Stage 1 first; Stage 2 adds independent dots from same latent line.\n\n")
    cat("Target slope:               ", round(attr(d2, "target_slope"), 4), "\n")
    cat("Stage 1 Deming slope:       ", round(attr(d1, "deming_slope"), 4), "\n")
    cat("Stage 2 Deming slope:       ", round(attr(d2, "deming_slope"), 4), "\n")
    cat("Stage 1 absolute error:     ", round(attr(d1, "stage1_error"), 4), "\n")
    cat("Stage 2 absolute error:     ", round(attr(d2, "stage2_error"), 4), "\n")
    cat("Stage 1 tries used:         ", attr(d2, "s1_try"), "\n")
    cat("Added-dot tries used:       ", attr(d2, "added_try"), "\n")
    cat("Stage 1 dots:               ", nrow(d1), "\n")
    cat("Stage 2 total dots:         ", nrow(d2), "\n")
    cat("New dots added at Stage 2:  ", attr(d2, "n_added"), "\n")
    cat("\nInterpretation:\n")
    cat("- Blue dashed line = latent target line, shown only for diagnosis.\n")
    cat("- Green line = Deming/orthogonal fit, used for regeneration diagnostics.\n")
    cat("- Stage 1 is not selected as a subset from Stage 2.\n")
    cat("- Stage 2 is accepted only when the full display's Deming slope stays close to target.\n")
  })
}

shinyApp(ui = ui, server = server)
