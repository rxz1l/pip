Rotate_Spatial_Data <- function(seurat_obj, 
                                 image_key = "sample1", 
                                 mode = "left", 
                                 scale_factor_type = "hires") {
  
  # 1. 提取图像对象
  vis_img <- seurat_obj@images[[image_key]]
  img_array <- vis_img@image
  H <- dim(img_array)[1] # 高度
  W <- dim(img_array)[2] # 宽度
  C <- dim(img_array)[3] # 通道
  
  # 2. 获取缩放因子并计算全分辨率参考维度
  sf <- vis_img@scale.factors[[scale_factor_type]]
  full_W <- W / sf
  full_H <- H / sf
  
  # 3. 准备坐标数据
  coords <- vis_img@coordinates
  old_x <- coords$imagecol
  old_y <- coords$imagerow
  
  # --- 核心变换逻辑 ---
  
  # 初始化新图像容器
  # 注意：90度旋转时宽高互换，180度或镜像时宽高不变
  if (mode %in% c("left", "right")) {
    new_img_array <- array(dim = c(W, H, C))
  } else {
    new_img_array <- array(dim = c(H, W, C))
  }

  for (i in 1:C) {
    channel <- img_array[, , i]
    
    if (mode == "left") {
      # 【左转 90°】: 转置 + 行逆序
      new_img_array[, , i] <- t(channel)[W:1, ]
      if (i == 1) { # 坐标只需计算一次
        coords$imagecol <- old_y
        coords$imagerow <- full_W - old_x
      }
      
    } else if (mode == "right") {
      # 【右转 90°】: 转置 + 列逆序
      new_img_array[, , i] <- t(channel)[, H:1]
      if (i == 1) {
        coords$imagecol <- full_H - old_y
        coords$imagerow <- old_x
      }
      
    } else if (mode == "180") {
      # 【旋转 180°】: 行逆序 + 列逆序
      new_img_array[, , i] <- channel[H:1, W:1]
      if (i == 1) {
        coords$imagecol <- full_W - old_x
        coords$imagerow <- full_H - old_y
      }
      
    } else if (mode == "mirror") {
      # 【水平镜像】: 列逆序
      new_img_array[, , i] <- channel[, W:1]
      if (i == 1) {
        coords$imagecol <- full_W - old_x
        coords$imagerow <- old_y
      }
    } else {
      stop("模式错误！请选择: 'left', 'right', '180', 或 'mirror'")
    }
  }
  
  # 4. 写回数据
  vis_img@image <- new_img_array
  vis_img@coordinates <- coords
  seurat_obj@images[[image_key]] <- vis_img
  
  message(paste("成功执行模式:", mode))
  return(seurat_obj)
}