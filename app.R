library(shiny)

# ---------- Deming regression ----------
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
                sqrt((Syy - lambda * Sxx)^2 + 4 * lambda * Sxy^2)) / (2 * Sxy)
  }
  
  intercept <- ybar - slope * xbar
  list(intercept = intercept, slope = slope)
}

# ---------- Generate one orthogonal-noise cloud ----------
generate_cloud <- function(target_slope, n_dots, noise_sd, bound, center_x = 0.5, center_y = 0.5) {
  t_vals <- seq(-0.45, 0.45, length.out = n_dots)
  
  v_line <- c(1, target_slope)
  v_line <- v_line / sqrt(sum(v_line^2))
  
  v_perp <- c(-target_slope, 1)
  v_perp <- v_perp / sqrt(sum(v_perp^2))
  
  d <- rnorm(n_dots, 0, noise_sd)
  tries <- 0
  while (any(abs(d) > bound) && tries < 5000) {
    d <- rnorm(n_dots, 0, noise_sd)
    tries <- tries + 1
  }
  if (tries >= 5000) return(NULL)
  
  x <- center_x + t_vals * v_line[1] + d * v_perp[1]
  y <- center_y + t_vals * v_line[2] + d * v_perp[2]
  
  data.frame(x = x, y = y, t = t_vals, d = d)
}

# ---------- More x-balanced subset sampling ----------
sample_balanced_subset <- function(data, n_stage1, n_bins = 3) {
  n <- nrow(data)
  
  if (n_stage1 >= n) return(seq_len(n))
  
  # Bin by x position
  breaks <- quantile(data$x, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
  # avoid duplicated breaks
  breaks <- unique(breaks)
  if (length(breaks) < 3) {
    return(sort(sample(seq_len(n), n_stage1, replace = FALSE)))
  }
  
  bins <- cut(data$x, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  bin_ids <- split(seq_len(n), bins)
  
  # baseline allocation across bins
  k <- length(bin_ids)
  alloc <- rep(floor(n_stage1 / k), k)
  remainder <- n_stage1 - sum(alloc)
  
  if (remainder > 0) {
    alloc[seq_len(remainder)] <- alloc[seq_len(remainder)] + 1
  }
  
  chosen <- integer(0)
  
  # first pass: sample up to allocation from each bin
  leftovers <- integer(0)
  for (i in seq_along(bin_ids)) {
    ids <- bin_ids[[i]]
    take <- min(length(ids), alloc[i])
    if (take > 0) {
      chosen <- c(chosen, sample(ids, take, replace = FALSE))
    }
    if (length(ids) > take) {
      leftovers <- c(leftovers, setdiff(ids, chosen))
    }
  }
  
  # fill remaining if needed
  still_needed <- n_stage1 - length(chosen)
  if (still_needed > 0) {
    pool <- setdiff(seq_len(n), chosen)
    chosen <- c(chosen, sample(pool, still_needed, replace = FALSE))
  }
  
  sort(chosen)
}

# ---------- Find Stage 2 + Stage 1 subset ----------
make_two_stage_trial <- function(
    target_slope = 0.10,
    n_stage2 = 20,
    n_stage1 = 8,
    noise_sd = 0.12,
    bound = 0.35,
    stage2_tol = 0.03,
    stage1_target_tol = 0.05,
    stage1_stage2_tol = 0.03,
    intercept_tol = 0.05,
    max_cloud_tries = 3000,
    max_subset_tries = 4000
) {
  
  if (n_stage1 >= n_stage2) stop("n_stage1 must be less than n_stage2.")
  
  for (cloud_try in 1:max_cloud_tries) {
    full <- generate_cloud(target_slope, n_stage2, noise_sd, bound)
    if (is.null(full)) next
    
    dem2 <- deming_fit(full$x, full$y)
    slope2 <- dem2$slope
    int2 <- dem2$intercept
    
    # Stage 2 must be close to target
    if (abs(slope2 - target_slope) > stage2_tol) next
    
    # Search for a Stage 1 subset within Stage 2
    for (subset_try in 1:max_subset_tries) {
      idx <- sample_balanced_subset(full, n_stage1, n_bins = 3)
      sub <- full[idx, , drop = FALSE]
      
      dem1 <- deming_fit(sub$x, sub$y)
      slope1 <- dem1$slope
      int1 <- dem1$intercept
      
      cond1 <- abs(slope1 - target_slope) <= stage1_target_tol
      cond2 <- abs(slope1 - slope2) <= stage1_stage2_tol
      cond3 <- abs(int1 - int2) <= intercept_tol
      
      if (cond1 && cond2 && cond3) {
        full$stage <- "Stage 2 only"
        full$included_stage1 <- FALSE
        full$included_stage1[idx] <- TRUE
        
        sub$stage <- "Stage 1"
        
        attr(full, "target_slope") <- target_slope
        attr(full, "dem2_slope") <- slope2
        attr(full, "dem2_intercept") <- int2
        attr(sub, "dem1_slope") <- slope1
        attr(sub, "dem1_intercept") <- int1
        attr(full, "cloud_try") <- cloud_try
        attr(full, "subset_try") <- subset_try
        
        return(list(stage1 = sub, stage2 = full))
      }
    }
  }
  
  NULL
}

# ---------- UI ----------
ui <- fluidPage(
  titlePanel("Two-Stage Deming-Stable Trial Generator"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("slope", "Target slope:",
                  min = -0.50, max = 0.50, value = 0.10, step = 0.01),
      sliderInput("noise", "Orthogonal noise SD:",
                  min = 0.01, max = 0.30, value = 0.10, step = 0.01),
      sliderInput("bound", "Noise bound:",
                  min = 0.05, max = 0.80, value = 0.30, step = 0.01),
      sliderInput("n2", "Stage 2 dots:",
                  min = 6, max = 40, value = 16, step = 1),
      sliderInput("n1", "Stage 1 dots:",
                  min = 3, max = 20, value = 7, step = 1),
      sliderInput("tol2", "Stage 2 slope tolerance:",
                  min = 0.005, max = 0.10, value = 0.03, step = 0.005),
      sliderInput("tol1target", "Stage 1-to-target tolerance:",
                  min = 0.005, max = 0.15, value = 0.05, step = 0.005),
      sliderInput("tol12", "Stage 1-to-Stage 2 tolerance:",
                  min = 0.005, max = 0.15, value = 0.03, step = 0.005),
      sliderInput("inttol", "Stage 1-to-Stage 2 intercept tolerance:",
                  min = 0.001, max = 0.20, value = 0.01, step = 0.001),
      numericInput("cloud_tries", "Max Stage 2 tries:", value = 3000, min = 100, step = 100),
      numericInput("subset_tries", "Max Stage 1 subset tries:", value = 4000, min = 100, step = 100),
      checkboxInput("show_target", "Show target line", TRUE),
      checkboxInput("show_deming1", "Show Stage 1 Deming line", TRUE),
      checkboxInput("show_deming2", "Show Stage 2 Deming line", TRUE),
      actionButton("regen", "Generate Trial")
    ),
    mainPanel(
      plotOutput("stage1plot", height = "320px"),
      plotOutput("stage2plot", height = "320px"),
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
      n_stage2 = input$n2,
      n_stage1 = input$n1,
      noise_sd = input$noise,
      bound = input$bound,
      stage2_tol = input$tol2,
      stage1_target_tol = input$tol1target,
      stage1_stage2_tol = input$tol12,
      intercept_tol = input$inttol,
      max_cloud_tries = input$cloud_tries,
      max_subset_tries = input$subset_tries
    )
  }, ignoreInit = FALSE)
  
  output$stage1plot <- renderPlot({
    tr <- trial()
    validate(need(!is.null(tr), "No valid trial found. Loosen tolerances, reduce noise, or increase tries."))
    
    d1 <- tr$stage1
    dem1s <- attr(d1, "dem1_slope")
    dem1i <- attr(d1, "dem1_intercept")
    
    plot(d1$x, d1$y,
         xlim = c(0, 1), ylim = c(0, 1),
         pch = 19, cex = 1.8,
         xlab = "", ylab = "",
         axes = FALSE, main = "Stage 1")
    box()
    
    if (input$show_target) {
      b <- input$slope
      a <- 0.5 - b * 0.5
      abline(a = a, b = b, col = "blue", lty = 2, lwd = 2.5)
    }
    
    if (input$show_deming1) {
      abline(a = dem1i, b = dem1s, col = "darkgreen", lwd = 2)
    }
  })
  
  output$stage2plot <- renderPlot({
    tr <- trial()
    validate(need(!is.null(tr), "No valid trial found. Loosen tolerances, reduce noise, or increase tries."))
    
    d2 <- tr$stage2
    dem2s <- attr(d2, "dem2_slope")
    dem2i <- attr(d2, "dem2_intercept")
    
    plot(d2$x, d2$y,
         xlim = c(0, 1), ylim = c(0, 1),
         pch = 19, cex = 1.6,
         xlab = "", ylab = "",
         axes = FALSE, main = "Stage 2")
    box()
    
    points(d2$x[!d2$included_stage1], d2$y[!d2$included_stage1], pch = 19, cex = 1.6, col = "gray75")
    points(d2$x[d2$included_stage1], d2$y[d2$included_stage1], pch = 19, cex = 1.8, col = "black")
    
    if (input$show_target) {
      b <- input$slope
      a <- 0.5 - b * 0.5
      abline(a = a, b = b, col = "blue", lty = 2, lwd = 2.5)
    }
    
    if (input$show_deming2) {
      abline(a = dem2i, b = dem2s, col = "darkgreen", lwd = 2)
    }
  })
  
  output$info <- renderPrint({
    tr <- trial()
    if (is.null(tr)) {
      cat("No valid trial found.\n")
      cat("Try loosening tolerances, lowering noise, raising bounds, or increasing tries.\n")
      return()
    }
    
    d1 <- tr$stage1
    d2 <- tr$stage2
    
    cat("Target slope:            ", round(attr(d2, "target_slope"), 4), "\n")
    cat("Stage 1 Deming slope:    ", round(attr(d1, "dem1_slope"), 4), "\n")
    cat("Stage 2 Deming slope:    ", round(attr(d2, "dem2_slope"), 4), "\n")
    cat("Slope difference:        ", round(abs(attr(d1, "dem1_slope") - attr(d2, "dem2_slope")), 4), "\n")
    cat("Stage 1 intercept:       ", round(attr(d1, "dem1_intercept"), 4), "\n")
    cat("Stage 2 intercept:       ", round(attr(d2, "dem2_intercept"), 4), "\n")
    cat("Intercept difference:    ", round(abs(attr(d1, "dem1_intercept") - attr(d2, "dem2_intercept")), 4), "\n")
    cat("Stage 2 generation try:  ", attr(d2, "cloud_try"), "\n")
    cat("Stage 1 subset try:      ", attr(d2, "subset_try"), "\n")
    cat("Stage 1 dots:            ", nrow(d1), "\n")
    cat("Stage 2 dots:            ", nrow(d2), "\n")
  })
}

shinyApp(ui = ui, server = server)