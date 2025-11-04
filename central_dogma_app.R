library(shiny)

# ----- Functions ------------------------------------------------------------

# Virtual gene: make a DNA string of length divisible by 3
gene_dna <- function(length, base_probs = c(0.25, 0.25, 0.25, 0.25)){
  if (length %% 3 != 0){
    stop("The 'length' must be divisible by 3 for codons.")
  }
  dna_vector <- sample(
    x = c("A", "T", "C", "G"),
    size = length,
    replace = TRUE,
    prob = base_probs
  )
  paste0(dna_vector, collapse = "")
}

# Transcription: DNA -> RNA (T -> U)
transcribe_dna <- function(dna){
  gsub(pattern = "T", replacement = "U", x = dna)
}

# Translation: RNA -> Protein (AA letters), stop at first STOP
translate_rna <- function(rna){
  # Genetic code (RNA codons)
  code <- c(
    UUU="F", UUC="F", UUA="L", UUG="L",
    CUU="L", CUC="L", CUA="L", CUG="L",
    AUU="I", AUC="I", AUA="I", AUG="M",
    GUU="V", GUC="V", GUA="V", GUG="V",
    UCU="S", UCC="S", UCA="S", UCG="S",
    CCU="P", CCC="P", CCA="P", CCG="P",
    ACU="T", ACC="T", ACA="T", ACG="T",
    GCU="A", GCC="A", GCA="A", GCG="A",
    UAU="Y", UAC="Y", UAA="*", UAG="*",
    CAU="H", CAC="H", CAA="Q", CAG="Q",
    AAU="N", AAC="N", AAA="K", AAG="K",
    GAU="D", GAC="D", GAA="E", GAG="E",
    UGU="C", UGC="C", UGA="*", UGG="W",
    CGU="R", CGC="R", CGA="R", CGG="R",
    AGU="S", AGC="S", AGA="R", AGG="R",
    GGU="G", GGC="G", GGA="G", GGG="G"
  )
  
  # Split into codons
  codons <- substring(rna, seq(1, nchar(rna), by = 3), seq(3, nchar(rna), by = 3))
  aas <- unname(code[codons])
  # Replace unknown/NA codons with "X"
  aas[is.na(aas)] <- "X"
  
  # Stop at first STOP (*)
  stop_idx <- match("*", aas, nomatch = 0)
  if (stop_idx > 0) aas <- aas[seq_len(stop_idx - 1)]
  
  paste0(aas, collapse = "")
}

# Base frequencies: counts and proportions
base_freqs <- function(seq, alphabet = c("A","C","G","T","U")){
  chars <- strsplit(seq, "")[[1]]
  bases <- intersect(alphabet, unique(chars))
  counts <- sapply(alphabet, function(b) sum(chars == b))
  props  <- if (length(chars) > 0) counts / length(chars) else counts
  data.frame(Base = alphabet, Count = counts, Proportion = round(props, 3))
}

# ----- UI -------------------------------------------------------------------

ui <- fluidPage(
  titlePanel("Virtual Central Dogma"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("n_bases", "Number of bases (multiple of 3):", min = 3, max = 300, value = 30, step = 3),
      numericInput("prob_A", "Probability of A", value = 0.25, min = 0, max = 1, step = 0.05),
      numericInput("prob_T", "Probability of T", value = 0.25, min = 0, max = 1, step = 0.05),
      numericInput("prob_C", "Probability of C", value = 0.25, min = 0, max = 1, step = 0.05),
      numericInput("prob_G", "Probability of G", value = 0.25, min = 0, max = 1, step = 0.05),
      actionButton("regen", "Generate new gene")
    ),
    mainPanel(
      h4("DNA"),
      verbatimTextOutput("dna"),
      h4("RNA (transcribed)"),
      verbatimTextOutput("rna"),
      h4("Protein (translated)"),
      verbatimTextOutput("protein"),
      h4("Base frequencies"),
      tableOutput("dna_freqs")
    )
  )
)

# ----- Server ---------------------------------------------------------------

server <- function(input, output, session){
  
  # Validate probabilities once when needed
  probs <- reactive({
    p <- c(input$prob_A, input$prob_T, input$prob_C, input$prob_G)
    if (any(is.na(p)) || any(p < 0) || any(p > 1))
      validate(need(FALSE, "Probabilities must be between 0 and 1."))
    s <- sum(p)
    validate(need(abs(s - 1) < 1e-8, "Probabilities must sum to 1."))
    p
  })
  
  # Generate DNA whenever the user clicks "Generate" or changes length/probs
  dna_seq <- reactive({
    # tie to inputs so changes also regenerate
    input$regen
    isolate({
      validate(need(input$n_bases %% 3 == 0, "Length must be divisible by 3."))
      gene_dna(length = input$n_bases, base_probs = probs())
    })
  })
  
  # Downstream reactive derivations
  rna_seq <- reactive( transcribe_dna(dna_seq()) )
  protein_seq <- reactive( translate_rna(rna_seq()) )
  
  # Outputs
  output$dna     <- renderText(dna_seq())
  output$rna     <- renderText(rna_seq())
  output$protein <- renderText({
    p <- protein_seq()
    if (nchar(p) == 0) "(Empty protein: translation stopped immediately or no codons)"
    else p
  })
  output$dna_freqs <- renderTable({
    base_freqs(dna_seq(), alphabet = c("A","C","G","T"))
  }, striped = TRUE, spacing = "s")
}

shinyApp(ui, server)
