# Load packages ----
# see here: https://community.rstudio.com/t/failing-to-deploy-shinyapp-depending-on-bioconductor-packages/6970/4
library(BiocManager)
options(repos = BiocManager::repositories())

library(shiny)
library(Biostrings)
library(magrittr)
library(DT)
library(janitor)
library(tidyverse)


# Source helper functions, load data -----
source("reverse_translation.R")



# User interface ----
ui <- fluidPage(
  titlePanel("Reverse translation tool"),
  
  sidebarLayout(
    sidebarPanel(
      "", 
      h5("1.) Paste your amino acid sequences in fasta format here:"), 
      textAreaInput(
        inputId = "sequence_input", 
        label = "", 
        placeholder = ">seq1\nMTRSDAAGH"
      ),
      "Alternatively: upload a fasta file: ", 
      fileInput("sequence_input_fasta", "", placeholder = ".fasta"), 
      
      actionButton(
        "submit", 
        "Reverse translate!"
      ), 
      
      fluidRow(
        h5("2.) Select codon usage table"), 
        column(3, 
          radioButtons("radio", label = "", 
                       choices = c("From species", "User defined"), 
                       selected = "From species")
        ), 
        column(9, 
          uiOutput("codon_choice")
        )
      ),
      
      fluidRow(
        h5("3.) Select restriction enzymes to exclude"), 
        checkboxGroupInput("enzymes", label = "",
                     choices = names(enzymes), 
                     selected = character(0), 
                     inline = TRUE),
      ), 
      
      fluidRow(
        h5("4.) Optional: Add additional DNA sequences you wish to exclude (one line per sequence)"), 
        textAreaInput(
          "exclusion_sites", label = "", value = NULL
        )
      ), 
      
    ), 
    mainPanel(
      h3("Reverse translation (codon optimized):"), 
      verbatimTextOutput("reverse_translation", placeholder = TRUE), 
      br(),
      downloadButton("download", "Download .fasta"), 
      br(), 
      h3("DNA sequences excluded in reverse translation:"), 
      verbatimTextOutput("cut_sites", placeholder = TRUE), 
      br(), 
      h3("Used codon usage table:"), 
      DTOutput("codon_usage_table")
    )
  ), 
  hr(), 
  div(
    class = "footer",
    includeHTML("footer.html")
  )
)

# Server logic ----
server <- function(input, output) {
  
  codon_usage_tab_active <- reactiveVal()
  codon_usage_tab_active_split <- reactive(
    split(codon_usage_tab_active(), codon_usage_tab_active()$aa_single)
  )
  
  retrieve_input_seqs <- eventReactive(
    input$submit, {
      if (!is.null(input$sequence_input_fasta$name)) {
        # check for file upload
        #browser()
        sequence_input <- readLines(con = input$sequence_input_fasta$datapath)
        message("Processing file input: ", sequence_input)
      } else {
        message("Processing text area input")
        sequence_input <- input$sequence_input
      }
      tmp_file <- tempfile(fileext = ".fasta")
      ## TO DO: remove spaces, points, digits
      writeLines(sequence_input, con = tmp_file)
      
      fasta_input <- readAAStringSet(tmp_file, format = "fasta")
      return(fasta_input)
    }
  )
  
  collect_cut_sites <- eventReactive(
    input$submit, {
      restriction_sites <- enzymes[input$enzymes]
      custom_sites <- str_split(input$exclusion_sites, pattern = "\n")[[1]]
      cut_sites <- c(restriction_sites, custom_sites)
      cut_sites <- cut_sites[nchar(cut_sites) != 0]
      if (length(cut_sites) == 0) {
        cut_sites <- NULL
      } else {
        names(cut_sites)[names(cut_sites) == ""] <- "(custom)"
      }
      return(cut_sites)
    }
  )
  
  generate_reverse_translation <- reactive({
    input_seqs <- retrieve_input_seqs()
    exclusion_sites <- collect_cut_sites()
    reverse_translated <- 
      lapply(input_seqs, 
             function(seq) {
               ## TO DO: passing names is not working atm
               reverse_translate(aa_seq = seq, aa_seq_name = names(seq), 
                                 codon_usage_tables = isolate(codon_usage_tab_active_split()),  
                                 exclude_sites = exclusion_sites)
             })
    reverse_translated <- DNAStringSet(reverse_translated)
    not_successful <- paste0(names(reverse_translated)[lengths(reverse_translated) == 0], collapse = "\n")
    if (nchar(not_successful) > 0) {
      msg <- paste0("Reverse translation under the provided constraints was not successful for the following sequences:\n", not_successful)
      warning(msg)
      showNotification(msg, duration = NULL, type = "warning")
    }
    reverse_translated <- as.character(reverse_translated)
    reverse_translated <- paste0(">", names(reverse_translated), "\n", reverse_translated, collapse = "\n")
    return(reverse_translated)
  })
  
  observeEvent(
    input$codon_choice, {
      message("observeEvent triggered by input$codon_choice")
      #browser()
      codon_usage_tab_active(codon_usage_tables[[input$codon_choice]])
  })
  
  observeEvent(
    input$codon_table_user, {
      message("observeEvent triggered by input$codon_table_user")
      codon_usage_tab_active(import_codon_usage_table(input$codon_table_user$datapath))
    })
  
  observeEvent(
    input$radio, {
      if (input$radio == "From species") {
        req(input$codon_choice)
        codon_usage_tab_active(codon_usage_tables[[input$codon_choice]])
      }
    }
  )
  
  output$reverse_translation <- renderText({
    generate_reverse_translation()
  })
  
  output$cut_sites <- renderPrint({
    my_out <- collect_cut_sites()
    if (is.null(my_out)) {
      return("(none)")
    } else {
      return(my_out)
    }
  })
  
  output$download <- downloadHandler(
    filename = "reverse_translation.fasta",
    content = function(file) {
      writeLines(generate_reverse_translation(), file)
    }
  )
  
  output$codon_choice <- renderUI({
    if (input$radio == "From species") {
      selectInput("codon_choice", label = '"rare codons excluded": codons with a frequency < 0.1 are only considered if sequence exclusion constraints cannot otherwise be satisfied',
                  #choices = c("S. cerevisiae", "E. coli", "H. sapiens"))
                  choices = list("S. cerevisiae" = "codon_usage_tbl_yeast", 
                                 "S. cerevisiae, rare codons excluded" = "codon_usage_tbl_yeast_rare_excluded", 
                                 "E. coli" = "codon_usage_tbl_ecoli", 
                                 "E. coli, rare codons excluded" = "codon_usage_tbl_ecoli_rare_excluded", 
                                 "H. sapiens" = "codon_usage_tbl_hsapiens", 
                                 "H. sapiens, rare codons excluded" = "codon_usage_tbl_hsapiens_rare_excluded"), 
                  selected = "codon_usage_tbl_yeast_rare_excluded")
    } else if (input$radio == "User defined") {
      fileInput(
        "codon_table_user", 
        label = "Provide a codon usage table in csv format with the following columns: AA_single_letter, AA_triple_letter, codon, codon_frequency", 
        accept = ".csv"
      )
    }
  })
  
  output$codon_usage_table <- 
    renderDT({
      codon_usage_tab_active()
    }, 
      options = list(pageLength = 10)
    )
  
}

# Run app ----
shinyApp(ui, server)
