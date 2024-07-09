library(bslib)
library(shiny)
library(bsicons)


options(shiny.maxRequestSize=5000*1024^2)

shinyApp(
  ui = page_navbar(
    title = icon("biohazard"),
    theme = bs_theme(),
    collapsible = TRUE,
    inverse = TRUE,
    footer = card(max_height="200px",
      class="bg-warning",
      card_header("Disclaimer"),
      card_body("IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR 
        ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
        TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
        OR THE USE OR OTHER DEALINGS IN THE SOFTWARE."
      )
    ),
    shinyjs::useShinyjs(),
    #tags$head(tags$style(".modal-dialog{ width:1000px}")),
    #tags$head(tags$style(".modal-body{ min-height:700px}")),
    
    nav_panel(title="VCF",
              inputPanel(
                title = "Generate a report from a VCF containing variants",
                selectInput("resistance_db","Resistance database",choices = c("HSV1","HSV2")),
                fileInput("vcf_file","VCF file",placeholder = "Upload a VCF file (with INFO fields DP,AF,BCSQ)"),
              ),
              actionButton("make_vcf_report","Generate Reports",icon = icon("cog")),
              downloadButton("vcf_report_html","HTML Report",icon = icon("download")),
              downloadButton("vcf_report_docx","DOCX Report",icon = icon("download"))
    ),
    nav_panel(title="BAM",
          p("Use this form to generate a report from a BAM file containing ONT reads aligned on a reference"),
          inputPanel(width="100%",
            fileInput("bam_file","BAM file",placeholder = "Upload BAM file of aligned reads on the reference",accept = ".bam",multiple = FALSE)
          ),
          shinyjs::disabled(actionButton("btn_view_bam_report","View BAM Report",icon = icon("view"),class = "btn-primary"))
          #shinyjs::disabled(downloadButton("btn_dl_bam_report_html","Download",icon = icon("download"),class="btn-primary"))
    )
    #tabPanel(title="FASTA",
               #p("Map given consensus FASTA file and detect variant relatively to the reference"),
               #selectInput("resistance_db","Resistance database",choices = c("HSV1","HSV2")),
               #textAreaInput("fasta_seq","DNA sequences (FASTA)",placeholder = ">UL23\nATGCTCATGACTGCATCTA...\n>UL30\nATGCTCATGACTGCATCTA...",height = "400px"),
               #actionButton("make_doc","Report.doc",icon = icon("download")),
               #shinycssloaders::withSpinner(textOutput("terminal_stdout",container = pre),type = 1)
    #)
  ),

  server = function(input, output, session) {
    
    bam_report_dir <- reactiveVal(NULL)
    
    observeEvent(input$bam_file,{
      validate(need(is.character(input$bam_file$datapath),"invalid BAM"))
      
      # Copy inputs to temp location
      tempdir <- tempdir()
      file.copy("notebooks/bam_report.Rmd", file.path(tempdir,"template.Rmd"), overwrite = TRUE)
      file.copy(input$bam_file$datapath, file.path(tempdir,input$bam_file$name), overwrite = TRUE)
      
      # Generate report
      withProgress(message="Generate report",{
        rmarkdown::render(
          file.path(tempdir,"template.Rmd"),
          output_file = file.path(tempdir,"bam_report.html"),
          params = list(input_bam_file=file.path(tempdir,input$bam_file$name)),
          envir = new.env(parent = globalenv())
        )
      })
      shinyjs::enable("btn_view_bam_report")
      shinyjs::enable("btn_dl_bam_report_html")
      bam_report_dir(tempdir)
    })
    
    observeEvent(input$btn_view_bam_report,{
      showModal(modalDialog(
        title = "BAM report",
        tags$iframe(srcdoc=readLines(file.path(bam_report_dir(),"bam_report.html")), style='width:100%;height:600px;'),
        easyClose = TRUE,size = "xl",
        footer = downloadButton("btn_dl_bam_report_html","Download",icon = icon("download"),class="btn-primary")
      ))
    })
    
    output$btn_dl_bam_report_html <- downloadHandler(
      filename = function(){paste0(input$bam_file$name,".html")},
      content = function(file) {file.copy(file.path(bam_report_dir(),"bam_report.html"),file)}
    )
    
    # output$terminal_stdout <- eventReactive(input$make_doc,{
    #   validate(need(is.character(input$vcf_file$datapath),"invalid fasta"))
    #   cmd <- sprintf("make %s 2>&1",input$vcf_file$datapath)
    #   print(cmd)
    #   withProgress(message="Generate report",{
    #     pipe(cmd) |> readLines() |> paste(collapse = "\n")
    #   })
    # })
  }
)
