// Test cases for perspective 3 point algorithms.
// Written by Simon Fuhrmann.

#include <gtest/gtest.h>

#include <vector>
#include <random>

#include "math/vector.h"
#include "math/matrix.h"
#include "math/matrix_qr.h"
#include "sfm/pose_p3p.h"

namespace
{
    void
    fill_colinear_points_and_directions (std::vector<math::Vec3d>* p,
        std::vector<math::Vec3d>* d)
    {
        math::Vec3d p1(-1.0, -1.0, 2.0);
        math::Vec3d p2(0.0, 0.0, 2.0);
        math::Vec3d p3(1.0, 1.0, 2.0);
        p->clear();
        p->push_back(p1);
        p->push_back(p2);
        p->push_back(p3);
        d->clear();
        d->push_back(p1.normalized());
        d->push_back(p2.normalized());
        d->push_back(p3.normalized());
    }

    void
    fill_test_points_and_directions (std::vector<math::Vec3d>* p,
        std::vector<math::Vec3d>* d)
    {
        math::Vec3d p1(-1.0, 1.0, 2.0);
        math::Vec3d p2(0.0, 0.0, 2.0);
        math::Vec3d p3(1.0, 1.0, 2.0);
        p->clear();
        p->push_back(p1);
        p->push_back(p2);
        p->push_back(p3);
        d->clear();
        d->push_back(p1.normalized());
        d->push_back(p2.normalized());
        d->push_back(p3.normalized());
    }

    void
    fill_groundtruth_data (std::vector<math::Vec3d>* p,
        std::vector<math::Vec3d>* d, math::Matrix<double, 3, 4>* s)
    {
        double mat[3 * 3];
        double mat_q[3 * 3];
        double mat_r[3 * 3];
        double rand_lo = -5.0, rand_hi = 5.0;
        
        srand (static_cast <unsigned> (time(0)));
        for (int i = 0; i < 9; ++i)
        {
            mat[i] = rand_lo + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(rand_hi-rand_lo)));
        }
    
        math::matrix_qr(mat, 3, 3, mat_q, mat_r, 1e-14);
        
        math::Matrix<double, 4, 4> pose(0.0);
        for(int i = 0; i < 3; ++i){
            for(int j = 0; j < 3; ++j){
                pose[i*4+j] = mat_q[i*3+j];
            }
        }
        
        for(int i = 0; i < 3; ++i){
            pose[i*4+3] = rand_lo + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(rand_hi-rand_lo)));
        }
        pose[15] = 1.0;
    
        math::Vec3d p1, p2, p3;
        
        for(int i = 0; i < 3; ++i){
            p1[i] = rand_lo + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(rand_hi-rand_lo)));
            p2[i] = rand_lo + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(rand_hi-rand_lo)));;
            p3[i] = rand_lo + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(rand_hi-rand_lo)));;
        }
        
        p->clear();
        p->push_back(p1);
        p->push_back(p2);
        p->push_back(p3);
        d->clear();
        d->push_back(pose.mult(p1, 1.0).normalized());
        d->push_back(pose.mult(p2, 1.0).normalized());
        d->push_back(pose.mult(p3, 1.0).normalized());
        *s = math::Matrix<double, 3, 4>(pose.begin());
    }
}

TEST(PoseP3PTest, NumSolutions)
{
    std::vector<math::Vec3d> points, directions;
    std::vector<math::Matrix<double, 3, 4> > solutions;
    fill_test_points_and_directions(&points, &directions);

    sfm::pose_p3p_kneip(points[0], points[1], points[2],
        directions[0], directions[1], directions[2], &solutions);
    EXPECT_EQ(4, solutions.size());

    fill_colinear_points_and_directions(&points, &directions);
    sfm::pose_p3p_kneip(points[0], points[1], points[2],
        directions[0], directions[1], directions[2], &solutions);
    EXPECT_EQ(0, solutions.size());
}

TEST(PoseP3PTest, GroundTruth1)
{
    std::vector<math::Vec3d> points, directions;
    math::Matrix<double, 3, 4> solution;
    std::vector<math::Matrix<double, 3, 4> > solutions;
    bool found_good_solutions = true;
    
    for(std::size_t t = 0; t < 1000000; ++t){
        fill_groundtruth_data(&points, &directions, &solution);

        sfm::pose_p3p_kneip(points[0], points[1], points[2],
            directions[0], directions[1], directions[2], &solutions);
        bool found_good_solution = false;
        for (std::size_t i = 0; i < solutions.size(); ++i)
        {
            if (solution.is_similar(solutions[i], 1e-10))
                found_good_solution = true;
        }
        found_good_solutions &= found_good_solution;
        if(!found_good_solutions) break;
    }
    
    EXPECT_TRUE(found_good_solutions);
}

