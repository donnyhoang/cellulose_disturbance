#make a quick visual representation of your hypothesis for presentations

library(ggplot2)


diversity <- c(10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
disturbance <- c(0.1, 0.2, 0.3, .4, .5, .6, .7, .8, .9, 1,0.1, 0.2, 0.3, .4, .5, .6, .7, .8, .9, 1)
substrate <- c("Cellulose","Cellulose","Cellulose","Cellulose","Cellulose","Cellulose","Cellulose","Cellulose","Cellulose","Cellulose",
               "Glucose","Glucose","Glucose","Glucose","Glucose","Glucose","Glucose","Glucose","Glucose","Glucose")
data <- data.frame(diversity, disturbance, substrate)

p1 <- ggplot(data, aes(x = disturbance, y = diversity, color = substrate)) +
  geom_line(lwd = 2) +
  scale_color_manual(values=c("#7fbf7b","#e9a3c9")) +
  theme_classic(base_size = 20) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
  guides(color=guide_legend("Substrate")) +
  labs(y = "Diversity", x = "Disturbance Frequency") +
  ggtitle("Predicted DDR")

p1

ggsave("hypoth2.png", device = "png", dpi = 700)

data2 <- subset(data, substrate == "Cellulose")
p2 <- ggplot(data2, aes(x = disturbance, y = diversity, color = substrate)) +
  geom_line(lwd = 2) +
  scale_color_manual(values=c("#7fbf7b")) +
  theme_classic(base_size = 20) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
  guides(color=guide_legend("Substrate")) +
  labs(y = "Diversity", x = "Disturbance Frequency") +
  ggtitle("Predicted DDR")

p2
ggsave("hypoth.png", device = "png", dpi = 700)
