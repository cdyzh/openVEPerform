library(shiny)
library(dplyr)
library(yogiroc)
library(rsconnect)

# Define UI for the application
ui <- fluidPage(
  titlePanel("VEPerform"),
  
  sidebarLayout(
    sidebarPanel(
      checkboxInput("upload_file", "Upload Your Own Reference Set", value = FALSE),
      conditionalPanel(
        condition = "input.upload_file == true",
        fileInput("file", "Upload CSV File", 
                  accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")), 
                  actionButton("infoButton", "Formatting Guide", class = "btn-info")
      ),
      selectizeInput("gene", "Select Gene Name:", choices = NULL, options = list(maxOptions = 1000)),
      checkboxInput("common_variant_filter", "Exclude Common Variants (gnomAD AF > 0.005)", value = TRUE),
      checkboxGroupInput("scores", "Select Scores to Include:",
                         choices = list(
                           "VARITY",
                           "REVEL",
                           "AlphaMissense"
                         ),
                         selected = c("VARITY", "REVEL", "AlphaMissense") # Check all by default
      ),
      actionButton("plotButton", "Generate PRC Plot"),
      downloadButton("downloadPlot", "Download PRC Plot"),
      helpText("Example genes: LDLR, TTN (no VARITY), ACADVL, DNM2, MYH7, SCN5A, etc."),
      helpText("If error occurs, it means there is not enough data to generate the PRC. Try using fewer predictors.")
    ),
    
    mainPanel(
      plotOutput("prcPlot", width = "600px", height = "600px"),
      textOutput("errorText")
    )
  )
)

server <- function(input, output, session) {
  plot_data <- reactiveVal(NULL)
  
  # Replace column names with standard names
  standard_colnames <- c(
    "base__hugo", 
    "gnomad__af", 
    "varity_r__varity_r", 
    "alphamissense__am_pathogenicity", 
    "revel__score", 
    "classification"
  )
  
  # Read inputted csv, else read existing csv
  prcdata <- reactive({
    if (input$upload_file) {
      req(input$file)
      df <- read.csv(input$file$datapath, stringsAsFactors = FALSE)
      # Automatically rename columns to match the standard column names
      colnames(df) <- standard_colnames
      df
    } else {
      read.table("preprocessed_full.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE) # 175166 rows after processing
    }
  })
  
  # Update the dropdown menu with gene names
  observe({
    df <- prcdata()
    if (!is.null(df)) {
      gene_names <- unique(df$base__hugo)
      updateSelectizeInput(session, "gene", choices = gene_names, selected = gene_names[1], server = TRUE)
    }
  })
  
  # Update which scores exist based on selected gene
  observe({
    req(input$gene)
    df <- prcdata()
    gene_data <- df %>% filter(base__hugo == input$gene)
    
    available_scores <- c()
    if (any(!is.na(gene_data$varity_r__varity_r))) {
      available_scores <- c(available_scores, "VARITY")
    }
    if (any(!is.na(gene_data$alphamissense__am_pathogenicity))) {
      available_scores <- c(available_scores, "AlphaMissense")
    }
    if (any(!is.na(gene_data$revel__score))) {
      available_scores <- c(available_scores, "REVEL")
    }
    
    updateCheckboxGroupInput(session, "scores", choices = available_scores, selected = available_scores)
  })
  
  # Plot
  observeEvent(input$plotButton, {
    df <- prcdata()
    gene_s <- input$gene
    exclude_common_variants <- input$common_variant_filter
    selected_scores <- input$scores
    
    names(df)[names(df) == "varity_r__varity_r"] <- "VARITY"
    names(df)[names(df) == "alphamissense__am_pathogenicity"] <- "AlphaMissense"
    names(df)[names(df) == "revel__score"] <- "REVEL"
    
    # Common variant filter
    if (exclude_common_variants) {
      df <- df[is.na(df$gnomad__af) | df$gnomad__af <= 0.005, ]
    }
    
    # Ordering is needed for yrobj to work
    df <- df[order(df$classification), ]
    
    # Gene filter
    prcfiltered <- df %>%
      filter(base__hugo == gene_s)
    
    B_org <- sum(prcfiltered$classification == TRUE & rowSums(!is.na(prcfiltered[selected_scores])) > 0)
    P_org <- sum(prcfiltered$classification == FALSE & rowSums(!is.na(prcfiltered[selected_scores])) > 0)
    
    # Generate the yrobj object
    yrobj <- yr2(truth = prcfiltered[["classification"]], scores = prcfiltered[selected_scores], high = rep(FALSE, length(selected_scores)))

    plot_data(list(
      yrobj = yrobj,
      lty_styles = c("dashed", "solid", "dashed")[1:length(selected_scores)],
      col_styles = c("purple", "cadetblue2", "orange")[1:length(selected_scores)],
      gene_s = gene_s,
      selected_scores = selected_scores,
      B_org = B_org,
      P_org = P_org
    ))
    
    # Attempt to generate plot, catch errors
    output$prcPlot <- renderPlot({
      plot_info <- plot_data()
      if (!is.null(plot_info)) {
        tryCatch({
          draw.prc(plot_info$yrobj, lty = plot_info$lty_styles, col = plot_info$col_styles, lwd = 2, balanced = TRUE, main = paste0(plot_info$gene_s, " PRCs for ", paste(plot_info$selected_scores, collapse = ", ")))
          abline(h = 90, lty = "dashed")
          legend("left", legend = c(paste("# of Pathogenic and Likely Pathogenic:", plot_info$P_org), paste("# of Benign and Likely Benign:", plot_info$B_org)), pch = 15, bty = "n")
        }, error = function(e) {
          showModal(modalDialog(
            title = 'Error',
            'Not enough data',
            easyClose = TRUE,
            footer = NULL
          ))
        })
      }
    }, width = 600, height = 600, res = 72)
  })
  
  # Guide on how to format user-inputted csv
  observeEvent(input$infoButton, {
    showModal(modalDialog(
      title = "Reference Set Format Information",
      HTML("Please ensure your own reference set is a CSV file and is formatted with columns in the following order:<br><br>
            <b>gene:</b> Gene names<br>
            <b>gnomad_af:</b> GnomAD allele frequency<br>
            <b>varity:</b> VARITY score<br>
            <b>alphamissense:</b> AlphaMissense score<br>
            <b>revel:</b> REVEL score<br>
            <b>classification:</b> Variant classification (TRUE for Pathogenic, FALSE for Benign)<br><br>
            Please make sure to keep all columns. You can put NA for gnomad if you do not need the common variant filter, and put NA for each predictor if you do not wish to include it. <br><br>
            <b>Example Table:</b><br>
            <table border='1'>
              <tr>
                <th>gene</th>
                <th>gnomad_af</th>
                <th>varity</th>
                <th>alphamissense</th>
                <th>revel</th>
                <th>classification</th>
              </tr>
              <tr>
                <td>LDLR</td>
                <td>NA</td>
                <td>0.85</td>
                <td>0.95</td>
                <td>NA</td>
                <td>TRUE</td>
              </tr>
              <tr>
                <td>LDLR</td>
                <td>NA</td>
                <td>0.65</td>
                <td>0.75</td>
                <td>NA</td>
                <td>FALSE</td>
              </tr>
            </table>"),
      easyClose = TRUE,
      footer = NULL
    ))
  })

  # Download plot
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("PRC_plot_", input$gene, ".png", sep = "")
    },
    content = function(file) {
      plot_info <- plot_data()
      if (!is.null(plot_info)) {
        png(file, width = 6, height = 6, units = "in", res = 72)
        draw.prc(plot_info$yrobj, lty = plot_info$lty_styles, col = plot_info$col_styles, lwd = 2, balanced = TRUE, main = paste0(plot_info$gene_s, " PRCs for ", paste(plot_info$selected_scores, collapse = ", ")))
        abline(h = 90, lty = "dashed")
        legend("left", legend = c(paste("# of Pathogenic and Likely Pathogenic:", plot_info$P_org), paste("# of Benign and Likely Benign:", plot_info$B_org)), pch = 15, bty = "n")
        dev.off()
      }
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
