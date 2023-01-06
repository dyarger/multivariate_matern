
library(ggplot2)
theme_set(theme_bw() + theme(text = element_text(size = 18)))
x <- seq(-5, 5, by = .02)
library(shiny)
library(fields)
library(dplyr)
source('../../code/multi_matern_source.R')
grid_info <- create_grid_info_1d(2^12, 25)
grid_info_2d <- create_grid_info_2d(2^10, 25)

ui <- fluidPage(
    tabsetPanel(tabPanel('Matern Covariance',
                         sidebarLayout(
                           sidebarPanel(
                             withMathJax(),
                             sliderInput("nu", label = shiny::withMathJax("$$\\nu$$"),
                                         min = .02,
                                         max = 5,
                                         value = .5),
                             sliderInput("a",
                                         label = shiny::withMathJax("$$a$$"),
                                         min = .02,
                                         max = 5,
                                         value = 1),
                             sliderInput("sigma2",
                                         shiny::withMathJax("$$\\sigma^2$$"),
                                         min = .02,
                                         max = 5,
                                         value = 1)
                           ),mainPanel(
                             plotOutput("distPlot")
                           ))),
                tabPanel('Cross-Covariance (d=1, real)',
                         sidebarLayout(
                           sidebarPanel(
                             sliderInput("nu1",
                                         shiny::withMathJax("$$\\nu_1$$"),
                                         min = .02,
                                         max = 5,
                                         value = .5),
                             sliderInput("nu2",
                                         shiny::withMathJax("$$\\nu_2$$"),
                                         min = .02,
                                         max = 5,
                                         value = .5),
                             sliderInput("a1",
                                         shiny::withMathJax("$$a_1$$"),
                                         min = .02,
                                         max = 5,
                                         value = 1),
                             sliderInput("a2",
                                         shiny::withMathJax("$$a_2$$"),
                                         min = .02,
                                         max = 5,
                                         value = 1),
                             sliderInput("sigma2_new",
                                         shiny::withMathJax("$$\\Sigma_{12}$$"),
                                         min = -5.5,
                                         max = 5.5,
                                         value = 1)
                           ),mainPanel(
                             plotOutput("distPlot2")
                           ))),
                tabPanel('Cross-Covariance (d=1, complex)',
                         sidebarLayout(
                           sidebarPanel(
                             sliderInput("nu12",
                                         shiny::withMathJax("$$\\nu_1$$"),
                                         min = .02,
                                         max = 5,
                                         value = .5),
                             sliderInput("nu22",
                                         shiny::withMathJax("$$\\nu_2$$"),
                                         min = .02,
                                         max = 5,
                                         value = .5),
                             sliderInput("a12",
                                         shiny::withMathJax("$$a_1$$"),
                                         min = .02,
                                         max = 5,
                                         value = 1),
                             sliderInput("a22",
                                         shiny::withMathJax("$$a_2$$"),
                                         min = .02,
                                         max = 5,
                                         value = 1),
                             sliderInput("sigma2_new2",
                                         shiny::withMathJax("$$\\Re(\\Sigma_{12})$$"),
                                         min = -5.5,
                                         max = 5.5,
                                         value = as.double(0)),
                             sliderInput("sigma2_im",
                                         shiny::withMathJax("$$\\Im(\\Sigma_{12})$$"),
                                         min = -5.5,
                                         max = 5.5,
                                         value = 1)
                           ),mainPanel(
                             plotOutput("distPlot3")
                           ))),
                # tabPanel('Cross-Covariance (d=2, real)',
                #          sidebarLayout(
                #            sidebarPanel(
                #              sliderInput("nu1_d2",
                #                          shiny::withMathJax("$$\\nu_1$$"),
                #                          min = .02,
                #                          max = 5,
                #                          value = .5),
                #              sliderInput("nu2_d2",
                #                          shiny::withMathJax("$$\\nu_2$$"),
                #                          min = .02,
                #                          max = 5,
                #                          value = .5),
                #              sliderInput("a1_d2",
                #                          shiny::withMathJax("$$a_1$$"),
                #                          min = .02,
                #                          max = 5,
                #                          value = 1),
                #              sliderInput("a2_d2",
                #                          shiny::withMathJax("$$a_2$$"),
                #                          min = .02,
                #                          max = 5,
                #                          value = 1),
                #              sliderInput("sigma2_d2",
                #                          shiny::withMathJax("$$\\Re(\\Sigma_{12})$$"),
                #                          min = -5.5,
                #                          max = 5.5,
                #                          value = as.double(1)),
                #              sliderInput("theta_star",
                #                          shiny::withMathJax("$$atan2(\\theta^*)$$"),
                #                          min = -round(pi, 3),
                #                          max = round(pi, 3),
                #                          value = as.double(0))
                #            ),mainPanel(
                #              plotOutput("d2")
                #            ))),
                tabPanel('Cross-Covariance (d=2, complex)',
                         sidebarLayout(
                           sidebarPanel(
                             sliderInput("nu1_d2_im",
                                         shiny::withMathJax("$$\\nu_1$$"),
                                         min = .02,
                                         max = 5,
                                         value = .5),
                             sliderInput("nu2_d2_im",
                                         shiny::withMathJax("$$\\nu_2$$"),
                                         min = .02,
                                         max = 5,
                                         value = .5),
                             sliderInput("a1_d2_im",
                                         shiny::withMathJax("$$a_1$$"),
                                         min = .02,
                                         max = 5,
                                         value = 1),
                             sliderInput("a2_d2_im",
                                         shiny::withMathJax("$$a_2$$"),
                                         min = .02,
                                         max = 5,
                                         value = 1),
                             sliderInput("sigma2_d2_im",
                                         shiny::withMathJax("$$\\Re(\\Sigma_{12})$$"),
                                         min = -5.5,
                                         max = 5.5,
                                         value = as.double(0)),
                             sliderInput("sigma2_d2_im_im",
                                         shiny::withMathJax("$$\\Im(\\Sigma_{12})$$"),
                                         min = -5.5,
                                         max = 5.5,
                                         value = as.double(1)),
                             sliderInput("theta_star_im",
                                         shiny::withMathJax("$$atan2(\\theta^*)$$"),
                                         min = -round(pi, 3),
                                         max = round(pi, 3),
                                         value = as.double(0)),
                             shiny::radioButtons(inputId = 'color_scale', 
                                                 label = 'Color Scale', 
                                                 choices = c('Bivariate' = 'gradient2',
                                                             'Rainbow'), selected = 'gradient2')
                           ),mainPanel(
                             plotOutput("d2_im")
                           )))),
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
      #print(length(x))
      #print(length(Matern(abs(x), alpha = input$a, smoothness = input$nu)))
      ggplot(data = data.frame('Lag'  =  x, `Covariance value`  = input$sigma2 * 
                                 Matern(abs(x), alpha = input$a, smoothness = input$nu)),
             aes(x = Lag, y = Covariance.value)) + 
        geom_line() + 
        labs(y = 'Covariance value')
    })
    
    output$distPlot2 <- renderPlot({
      # generate bins based on input$bins from ui.R
      df <- fft_1d(nu1 = input$nu1, nu2 = input$nu2, a1 = input$a1,
                   a2 = input$a2, re = input$sigma2_new, im = 0,  grid_info = grid_info)
      x_res <- approx(x = df[,1], y = df[,2], xout = x)$y
      ggplot(data = data.frame('Lag'  =  x, `Covariance value`  = x_res),
             aes(x = Lag, y = Covariance.value)) + 
        geom_line() + 
        labs(y = 'Cross-covariance value')
    })
    
    output$distPlot3 <- renderPlot({
      # generate bins based on input$bins from ui.R
      df <- fft_1d(nu1 = input$nu12, nu2 = input$nu22, a1 = input$a12,
                   a2 = input$a22, re = input$sigma2_new2, im = input$sigma2_im,  grid_info = grid_info)
      x_res <- approx(x = df[,1], y = df[,2], xout = x)$y
      ggplot(data = data.frame('Lag'  =  x, `Covariance value`  = x_res),
             aes(x = Lag, y = Covariance.value)) + 
               geom_line() + 
        labs(y = 'Cross-covariance value')
    })
    
    output$d2_im <- renderPlot({
      Psi_fun_im <- function(theta) {
        if (input$theta_star_im < 0) {
          ifelse(theta > input$theta_star_im & theta < input$theta_star_im + pi, 1, -1)
        } else {
          ifelse(theta > input$theta_star_im | theta < input$theta_star_im - pi, 1, -1)
        }
      }
      if (input$color_scale == 'gradient2') {
        color_scale <- scale_fill_gradient2() 
      } else {
        color_scale <- scale_fill_gradientn(colors = rev(rainbow(10))) 
      }
      Delta <- function(x) complex(real = input$sigma2_d2_im,
                                   imaginary = input$sigma2_d2_im_im*Psi_fun_im(x))
      # generate bins based on input$bins from ui.R
      df <- fft_2d(nu1 = input$nu1_d2_im, nu2 = input$nu2_d2_im, a1 = input$a1_d2_im,
                   a2 = input$a2_d2_im, Delta = Delta, 
                   Psi = Psi_fun_im, grid_info = grid_info_2d) %>%
        filter(abs(Var1) < 5, abs(Var2) < 5)
      ggplot(data = df,
             aes(x = Var1, y = Var2, fill = val, color = val)) + 
        geom_raster() + 
        coord_equal() + 
        color_scale + 
        labs(color = 'Cross-covariance value', fill = 'Cross-covariance value')
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
