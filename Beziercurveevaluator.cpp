#include "BezierCurveEvaluator.h"
#include <assert.h>
#include "mat.h"
#include "vec.h"

#define NUM_SEG 20

void BezierCurveEvaluator::evaluateCurve(const std::vector<Point>& ptvCtrlPts,
	std::vector<Point>& ptvEvaluatedCurvePts,
	const float& fAniLength,
	const bool& bWrap) const
{
	ptvEvaluatedCurvePts.clear();
	const Mat4d ker(
		-1, 3, -3, 1,
		3, -6, 3, 0,
		-3, 3, 0, 0,
		1, 0, 0, 0);
	std::vector<Point> ptvCtrlPtsCpy(ptvCtrlPts);
	if (bWrap)
	{
		ptvCtrlPtsCpy.push_back(Point(ptvCtrlPts.front().x + fAniLength, ptvCtrlPts.front().y));
		ptvCtrlPtsCpy.insert(ptvCtrlPtsCpy.begin(), Point(ptvCtrlPts.back().x - fAniLength, ptvCtrlPts.back().y));
	}

	int iter;

	bool is_wrapped = false;

	for (iter = 0; iter < ((int)ptvCtrlPtsCpy.size() - 3); iter += 3)
	{
		Vec4d ctrl_x(ptvCtrlPtsCpy[iter].x, ptvCtrlPtsCpy[iter + 1].x,
			ptvCtrlPtsCpy[iter + 2].x, ptvCtrlPtsCpy[iter + 3].x);
		Vec4d ctrl_y(ptvCtrlPtsCpy[iter].y, ptvCtrlPtsCpy[iter + 1].y,
			ptvCtrlPtsCpy[iter + 2].y, ptvCtrlPtsCpy[iter + 3].y);
		for (int i = 0; i < NUM_SEG; ++i)
		{
			double t = i / (double)NUM_SEG;
			Vec4d para(t*t*t, t*t, t, 1);
			Point eval_pt((float)(para*ker*ctrl_x), (float)(para*ker*ctrl_y));
			if (eval_pt.x > fAniLength && bWrap)
			{
				float x_mod = eval_pt.x - fAniLength;
				if (!is_wrapped)
				{
					Point prev_pt(ptvEvaluatedCurvePts.back());
					float mid_pt_y = (eval_pt.y*(fAniLength-prev_pt.x)+ prev_pt.y*(eval_pt.x-fAniLength))
						/(eval_pt.x - prev_pt.x);
					ptvEvaluatedCurvePts.push_back(Point(0, mid_pt_y));
					ptvEvaluatedCurvePts.push_back(Point(fAniLength, mid_pt_y));
					is_wrapped = true;
				}
				eval_pt.x = x_mod;
			}
			ptvEvaluatedCurvePts.push_back(eval_pt);
		}
		ptvEvaluatedCurvePts.push_back(ptvCtrlPtsCpy[iter + 3]);
	}

	for (; iter < (int)ptvCtrlPts.size(); iter++) ptvEvaluatedCurvePts.push_back(ptvCtrlPts[iter]);
	if (!bWrap)
	{
		ptvEvaluatedCurvePts.push_back(Point(0, ptvCtrlPts.front().y));
		ptvEvaluatedCurvePts.push_back(Point(fAniLength, ptvCtrlPts.back().y));
	}
	else if (!is_wrapped)
	{
		float mid_y = (ptvCtrlPts.back().y*ptvCtrlPts.front().x
			+ ptvCtrlPts.front().y*(fAniLength - ptvCtrlPts.back().x))
			/ (ptvCtrlPts.front().x + fAniLength - ptvCtrlPts.back().x);
		ptvEvaluatedCurvePts.push_back(Point(0, mid_y));
		ptvEvaluatedCurvePts.push_back(Point(fAniLength, mid_y));
	}
}