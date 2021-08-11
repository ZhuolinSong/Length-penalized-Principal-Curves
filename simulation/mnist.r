library(mnist)
data(mnist)
str(mnist)
   
train <- mnist$train$x
label <- mnist$train$y
rotate <- function(mat) {
    t(mat[nrow(mat):1, , drop = FALSE])
}

get_image <- function(i) {
    image_1 <- which(train[i, ] != 0)
    img_x <- image_1 %% 28
    img_x[which(img_x == 0)] <- 28
    img_y <- -floor(image_1 / 28)

    img_x <- img_x / 28
    img_y <- img_y / 28 + 1

    out <- list(x = img_x,
                y = img_y)
    return(out)
}

par(mfrow = c(3, 3))
for (i in seq(9)) {
    print(label[i])
    plot(get_image(i))
}


dev.new()
image(rotate(t(matrix((train[3, ]), nrow = 28, byrow = F))))
