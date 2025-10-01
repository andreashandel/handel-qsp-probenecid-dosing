library(flowdiagramr)

variables = c("U","I","V","A","F","S")
flows = list(U_flows = c("-b*U*V"),
             I_flows = c("b*U*V","-cI*I","-kA*A*I"),
             V_flows = c("(1-fV)*p*I/(1+kF*F)", "-cV*V"),
             A_flows = c("F*V/(F*V+hF)","gA*A"),
             F_flows = c("(1-fF)*gF*V/(V+hV)*(fmax-F)","-cF*F"),
             S_flows = c("gs*F","-cS*S")
)
model <- list(variables = variables, flows = flows)
layout = list(varlocations = matrix(c("U","I","V",
                                      "A","F","S"),
                                    nrow = 2, byrow = TRUE),
              varspace_x_size = 2.5
)


# make list structure for the diagram with defaults
# this won't look good, so we update below
dlist <- flowdiagramr::prepare_diagram(model,layout)

# move arrows and text around
update_settings = list(flow_label_text = c(i_bUV = "b*U*V",
                                           g_FVFVhF = "F*V/(F*V+hF)",
                                           i_1fVpI1kFF = "(1-fV)*p*I/(1+kF*F)"),
  flow_xstart = c(i_FVFVhF = 1.5,
                  e_kAAI    = 2,
                  g_FVFVhF = -3.5),
  flow_xend = c(i_FVFVhF = 0.5,
                i_kAAI    = -1),
  flow_ystart = c(e_cII = 1,
                  i_FVFVhF = 0.5,
                  g_FVFVhF = -2),
  flow_yend = c(e_cII = 1.5,
                i_FVFVhF = 0.5,
                i_kAAI    = 0.5),
  flow_curvature = c(i_kAAI = 0,
                     i_1fVpI1kFF = 0),
  flow_xlabel = c( e_cII = 0,
                   i_FVFVhF = 0.1,
                   i_1fFgFVVhVfmaxF = 0.3,
                   i_1fVpI1kFF = 0.4,
                   i_bUV = - 2.7,
                   i_kAAI    = - 0.5,
                   g_FVFVhF = - 1.8),
  flow_ylabel = c( e_cII = 1.7,
                   i_FVFVhF  = 0.5,
                   i_1fFgFVVhVfmaxF = -0.6,
                   i_1fVpI1kFF = 1.3,
                   i_bUV = - 1.5,
                   i_kAAI    = -0.1,
                   g_gsF = 0.1,
                   g_FVFVhF = - 0.9),
  flow_show_arrow = c(i_FVFVhF = FALSE,
                      i_kAAI = TRUE,
                      i_1fFgFVVhVfmaxF = TRUE,
                      m_1fFgFVVhVfmaxF = TRUE,
                      g_FVFVhF   = TRUE,
                      e_kAAI    = FALSE,
                      i_1fVpI1kFF = TRUE)
)


# update diagram with updated settings
dlist2 <- flowdiagramr::update_diagram(dlist,diagram_settings = update_settings)
# generate the updated diagram
diag3 <- flowdiagramr::make_diagram(dlist2)
plot(diag3)
# write the current diagram code to an R script file, overwrite if it exists
#flowdiagramr::write_diagram(diagram_list = dlist2, filename = 'model_diagram_code.R', always_overwrite = TRUE)
ggplot2::ggsave(here::here("assets","model-diagram.png"),diag3, width = 8, height = 5, units = "in", dpi = 300)



