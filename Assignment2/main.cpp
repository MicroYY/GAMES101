// clang-format off
#include <iostream>
#include <opencv2/opencv.hpp>
#include "rasterizer.hpp"
#include "global.hpp"
#include "Triangle.hpp"

constexpr double MY_PI = 3.1415926;
#define ANGLE2RADIAN(angle) angle * MY_PI / 180

Eigen::Matrix4f get_view_matrix(Eigen::Vector3f eye_pos)
{
    Eigen::Matrix4f view = Eigen::Matrix4f::Identity();

    Eigen::Matrix4f translate;
    translate << 1,0,0,-eye_pos[0],
                 0,1,0,-eye_pos[1],
                 0,0,1,-eye_pos[2],
                 0,0,0,1;

    view = translate*view;

    return view;
}

Eigen::Matrix4f get_model_matrix(float angle_x, float angle_y, float angle_z)
{
    //Eigen::Matrix4f model = Eigen::Matrix4f::Identity();

    Eigen::Matrix4f rotateMatrix_x = Eigen::Matrix4f::Identity();
    rotateMatrix_x << 1, 0, 0, 0,
        0, cos(ANGLE2RADIAN(angle_x)), -sin(ANGLE2RADIAN(angle_x)), 0,
        0, sin(ANGLE2RADIAN(angle_x)), cos(ANGLE2RADIAN(angle_x)), 0,
        0, 0, 0, 1;

    Eigen::Matrix4f rotateMatrix_y = Eigen::Matrix4f::Identity();
    rotateMatrix_y << cos(ANGLE2RADIAN(angle_y)), 0, sin(ANGLE2RADIAN(angle_y)), 0,
        0, 1, 0, 0,
        -sin(ANGLE2RADIAN(angle_y)), 0, cos(ANGLE2RADIAN(angle_y)), 0,
        0, 0, 0, 1;

    Eigen::Matrix4f rotateMatrix_z = Eigen::Matrix4f::Identity();
    rotateMatrix_z << cos(ANGLE2RADIAN(angle_z)), -sin(ANGLE2RADIAN(angle_z)), 0, 0,
        sin(ANGLE2RADIAN(angle_z)), cos(ANGLE2RADIAN(angle_z)), 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1;

    return rotateMatrix_x * rotateMatrix_y * rotateMatrix_z;
}

Eigen::Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio, float zNear, float zFar)
{
    // TODO: Copy-paste your implementation from the previous assignment.
    Eigen::Matrix4f projection;

    Eigen::Matrix4f orthographic = Eigen::Matrix4f::Identity();

    float t = tan(ANGLE2RADIAN(eye_fov / 2)) * abs(zNear);
    float r = aspect_ratio * t;
    float l = -r;
    float b = -t;
    float n = zNear;
    float f = zFar;

    orthographic << 2 / (r - l), 0, 0, -(r + l) / 2,
        0, 2 / (t - b), 0, -(t + b) / 2,
        0, 0, 2 / (n - f), -(n + f) / 2,
        0, 0, 0, 1;

    Eigen::Matrix4f persp2ortho = Eigen::Matrix4f::Identity();
    persp2ortho << n, 0, 0, 0,
        0, n, 0, 0,
        0, 0, n + f, -n * f,
        0, 0, 1, 0;

    projection = orthographic * persp2ortho;

    return projection;
}

int main(int argc, const char** argv)
{
    float angle_x = 0;
    float angle_y = 0;
    float angle_z = 0;
    bool command_line = false;
    std::string filename = "output.png";

    if (argc == 2)
    {
        command_line = true;
        filename = std::string(argv[1]);
    }

    rst::rasterizer r(700, 700);

    Eigen::Vector3f eye_pos = {0,0,5};


    std::vector<Eigen::Vector3f> pos
            {
                    {2, 0, -2},
                    {0, 2, -2},
                    {-2, 0, -2},
                    {3.5, -1, -5},
                    {2.5, 1.5, -5},
                    {-1, 0.5, -5}
            };

    std::vector<Eigen::Vector3i> ind
            {
                    {0, 1, 2},
                    {3, 4, 5}
            };

    std::vector<Eigen::Vector3f> cols
            {
                    {217.0, 238.0, 185.0},
                    {217.0, 238.0, 185.0},
                    {217.0, 238.0, 185.0},
                    {185.0, 217.0, 238.0},
                    {185.0, 217.0, 238.0},
                    {185.0, 217.0, 238.0}
            };

    auto pos_id = r.load_positions(pos);
    auto ind_id = r.load_indices(ind);
    auto col_id = r.load_colors(cols);

    int key = 0;
    int frame_count = 0;

    if (command_line)
    {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle_x, angle_y, angle_z));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, col_id, rst::Primitive::Triangle);
        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::cvtColor(image, image, cv::COLOR_RGB2BGR);

        cv::imwrite(filename, image);

        return 0;
    }

    while(key != 27)
    {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle_x, angle_y, angle_z));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, col_id, rst::Primitive::Triangle);

        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::cvtColor(image, image, cv::COLOR_RGB2BGR);
        cv::imshow("image", image);
        key = cv::waitKey(10);

        std::cout << "frame count: " << frame_count++ << '\n';
    }

    return 0;
}
// clang-format on