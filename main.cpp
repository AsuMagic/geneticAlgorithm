#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <ctime>
#include <random>
#include <onidev.h>

const od::Color sky_lightblue(203, 219, 252);
const od::Color sky_blue(120, 160, 255);

inline double dectorad(double x)
{
	return (x + 90) * 3.141 / 180.;
}

class target
{
public:
	target(od::Vec2 position) : position(position) {}

	void draw()
	{
		od::drawSpriteScaled("target", 0, position.x, position.y, 3.);
	}

	od::Vec2 getPosition() const { return position; }

private:
	od::Vec2 position;
};

class canonball
{
public:
	canonball(od::Vec2 position, float force, float angle, float zoomRatio, od::Vec2 targetPos, od::Vec2 canonPos) : position(position),
																													 forcex(sin(dectorad(angle)) * force * zoomRatio), forcey(cos(dectorad(angle)) * force * zoomRatio),
																													 angle(angle),
																													 zoomRatio(zoomRatio),
																													 distCanonTarget(od::pointDistance(targetPos.x, targetPos.y, canonPos.x, canonPos.y)),
																													 spriteid(od::spriteIndex("canonball")),
																													 canonPos(canonPos) {}

	void tick(float deltaTime)
	{
		if (!hasReachedLimits && (distCanonTarget > od::pointDistance(position.x, position.y, canonPos.x, canonPos.y)))
		{
			forcey += 0.3 * deltaTime;
			spriteAngle -= 1.5f * od::pointDistance(0.f, 0.f, forcex, forcey) * deltaTime;

			position.x += forcex * deltaTime;
			position.y += forcey * deltaTime;
		}
		else
		{
			hasReachedLimits = true;
		}
	}

	void draw(float alpha)
	{
		glColor4f(1.f, 1.f, 1.f, alpha);
		od::drawSpriteExt(spriteid, 0, position.x, position.y, 3., 3., spriteAngle);
		glColor4ub(255, 255, 255, 255);
	}

	bool objectiveReached() { return hasReachedLimits; }

	od::Vec2 getPosition() const { return position; }

private:
	od::Vec2 position;
	float forcex, forcey, angle, zoomRatio, spriteAngle = 0.f, distCanonTarget;
	int spriteid;
	od::Vec2 canonPos;
	bool hasReachedLimits = false;
};

// @TODO : Templates?
struct gene
{
	float data, minData, maxData, variationRange;
};

struct DNAObject
{
	float score = 0;
	std::vector<gene> genes;
};

bool scoreSorting(DNAObject obj1, DNAObject obj2) { return obj1.score > obj2.score; }

class geneticAlgorithm
{
public:
	geneticAlgorithm() : mt(time(NULL)) {}

	virtual void newGeneration()
	{
		currentScore = 0.;
		for (int i = 0; i < std::min(10, static_cast<int>(objects.size())); ++i)
			currentScore += objects[i].score / objects.size();

		lastScore = 0.;
		for (int i = 0; i < std::min(10, static_cast<int>(lastGenObjects.size())); ++i)
			lastScore += lastGenObjects[i].score / lastGenObjects.size();

		if (lastScore > currentScore)
		{
			//objects = lastGenObjects;
		}
		else
		{
			lastGenObjects = objects;
			++gen;
		}

		std::sort(objects.begin(), objects.end(), scoreSorting);
		objects.resize(std::min(static_cast<unsigned int>(5), objects.size()));

		for (gene& i : basegenes)
		{
			i.variationRange *= 0.95;
		}

		for (unsigned int i = 0; i < objects.size() - 1; i += 2)
		{
			DNAObject children;
			children.genes.resize(objects[i].genes.size());
			for (unsigned int j = 0; j < objects[i].genes.size(); ++j)
			{
				children.genes[j].data = (objects[i].genes[j].data + objects[i+1].genes[j].data) / 2.f;
			}

			objects.push_back(children);
		}

		for (unsigned int i = 0; i < objects.size(); i++)
		{
			for (unsigned int j = 0; j < objects[i].genes.size(); ++j)
			{
				std::uniform_real_distribution<float> dist(-basegenes[j].variationRange, basegenes[j].variationRange);
				float randval = dist(mt);
				objects[i].genes[j].data = std::min(basegenes[j].maxData, std::max(objects[i].genes[j].data + randval, basegenes[j].minData));
			}
		}
	}

	void setScoreFor(int index, float score) { objects[index].score = score; }
	float getScoreFor(int index) const { return objects[index].score; }

	double getLastIterationScore() { return currentScore; }
	double getLastGenerationScore() { return lastScore; }

	virtual int addDNAObject(DNAObject obj)
	{
        objects.push_back(obj);
        return objects.size() - 1;
	}

	int addRandomDNAObject()
	{
		DNAObject obj;
		obj.genes.resize(basegenes.size());
		for (unsigned int i = 0; i < basegenes.size(); ++i)
		{
			std::uniform_real_distribution<float> dist(basegenes[i].minData, basegenes[i].maxData);
            obj.genes[i].data = dist(mt);
		}

		return addDNAObject(obj);
	}

	const DNAObject& getObjectData(int index)
	{
		return objects[index];
	}

	unsigned int getGeneration() const { return gen; }
	unsigned int getObjectCount() const { return objects.size(); }

protected:
	unsigned int gen = 0;
	std::vector<gene> basegenes;
	std::vector<DNAObject> objects, lastGenObjects;

	double currentScore = 0., lastScore = 0.;
	std::mt19937 mt;
};

class canonGA : public virtual geneticAlgorithm
{
public:
	canonGA(od::Vec2 canonPos, od::Vec2 targetPos, float zoomRatio) : canonPos(canonPos), targetPos(targetPos), zoomRatio(zoomRatio)
	{
		basegenes.push_back({20.f, 7.5f, 20.f, 0.2f}); // fireIntensity
		basegenes.push_back({45.f, 25.f, 70.f, 1.f}); // angle
	}

	void newGeneration() override
	{
		_balls.clear();

		_scoreboard.resize(objects.size());
		geneticAlgorithm::newGeneration();

		for (unsigned int i = 0; i < objects.size(); ++i)
		{
			_scoreboard[i] = objects[i];
			_balls.push_back(canonball(canonPos, objects[i].genes[force].data, objects[i].genes[angle].data, zoomRatio, targetPos, canonPos));
		}

		_scoreevol.push_back(currentScore);
	}

	int addDNAObject(DNAObject obj) override
	{
		_balls.push_back(canonball(canonPos, obj.genes[force].data, obj.genes[angle].data, zoomRatio, targetPos, canonPos));
		return geneticAlgorithm::addDNAObject(obj);
	}

	void tick()
	{
		bool hasActive = false;

		for (unsigned int i = 0; i < balls.size(); ++i)
		{
			_balls[i].tick(od::deltaTime());

			float alpha = 60.f;
			if (_balls[i].objectiveReached())
			{
				setScoreFor(i, 20000 - od::pointDistance(targetPos.x, targetPos.y, _balls[i].getPosition().x, _balls[i].getPosition().y));
				alpha = 1.f;
			}
			else
			{
				hasActive = true;
			}

			_balls[i].draw(alpha / _balls.size());
		}

		if (!hasActive)
			newGeneration();
	}

	enum canonGenes
	{
		force = 0,
		angle
	};

	const std::vector<DNAObject>& scoreboard = _scoreboard;
	const std::vector<canonball>& balls = _balls;
	const std::vector<float>& scoreevol = _scoreevol;

private:
	od::Vec2 canonPos, targetPos;

	std::vector<canonball> _balls;
	std::vector<DNAObject> _scoreboard;
	std::vector<float> _scoreevol;

	float zoomRatio;
};

// Workaround around MinGW-related g++ bug #52015 bugzilla entry
#ifdef __MINGW32__
template<typename T>
std::string to_string(T type) // @TODO : Probably don't use string streams because they're heavy. To remove when the bug gets fixed in a stable version.
{
	std::stringstream ss;
	ss << type;
	return ss.str();
}
#endif

// @TODO : Way too much verbose, shorten it somehow.
void drawRectangleColor(od::Vec2 p1, od::Vec2 p2, od::Color c1, od::Color c2, od::Color c3, od::Color c4)
{
	glBegin(GL_TRIANGLE_STRIP);

    glColor4ub(c1.r, c1.g, c1.b, c1.a); glVertex2d(p1.x, p1.y); // top left
    glColor4ub(c2.r, c2.g, c2.b, c2.a); glVertex2d(p2.x, p1.y); // top right
    glColor4ub(c3.r, c3.g, c3.b, c3.a); glVertex2d(p2.x, p2.y); // bottom right
    glColor4ub(c1.r, c1.g, c1.b, c1.a); glVertex2d(p1.x, p1.y); // top left
    glColor4ub(c4.r, c4.g, c4.b, c4.a); glVertex2d(p1.x, p2.y); // bottom left

    glColor4ub(255, 255, 255, 255);

    glEnd();
}

void drawTextShadowed(od::Font& fnt, float x, float y, const std::string& str, float offsetx = 1, float offsety = 1, int halign = 0, int valign = 0, od::Color brightcolor = od::Color(240, 240, 240, 255), od::Color shadowcolor = od::Color(0, 0, 0, 120))
{
	// Shadow
	glColor4ub(shadowcolor.r, shadowcolor.g, shadowcolor.b, shadowcolor.a);
	fnt.draw(x + offsetx, y + offsety, str, halign, valign);

	// Regular
	glColor4ub(brightcolor.r, brightcolor.g, brightcolor.b, brightcolor.a);
	fnt.draw(x, y, str, halign, valign);
}

int main()
{
	using namespace od;
	using namespace std;

	Window win(displayGetWidth(), displayGetHeight(), "Canon genetic algorithm");
	win.setStyle(win.Style::Surface);
	win.setSynchronization(true);
    framerateSetDelta(60);

    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

    assetsLoadSprites("assets.xml");
    int canonSpriteID = spriteIndex("canon");

	float canonx = canonSpriteID / 2 + 96, canony = win.height() - 24;
    float zoomRatio = (displayGetWidth() / 800) * (displayGetHeight() / 600);

    Font fnt("roboto.ttf", 18 * zoomRatio, true);

    mt19937 mt(time(NULL));
    uniform_int_distribution<int> randx(win.width() / 2, win.width() - 64);
    uniform_int_distribution<int> randy(64, win.height() - 64);

	target t(Vec2(randx(mt), randy(mt)));

    canonGA algorithm(Vec2(canonx, canony), t.getPosition(), zoomRatio);

	for (int i = 0; i < 10; ++i)
		algorithm.addRandomDNAObject();

    while(win.open())
    {
        win.ioHandle();

        if(keyPressed(vk_escape))
            break;

        drawClear(Color(0, 0, 0).toRgb());

        win.updateView();
        drawRectangleColor(Vec2::Zero, Vec2(win.width(), win.height()), sky_lightblue, sky_lightblue, sky_blue, sky_blue);

		drawSpriteScaled(canonSpriteID, 1, canonx, canony, 3.);
		drawSpriteExt(canonSpriteID, 2, canonx, canony, 3., 3., 25);
		drawSpriteScaled(canonSpriteID, 0, canonx, canony, 3.);

		algorithm.tick();

		t.draw();

		int textOffset = 8 * zoomRatio;

		string geneScoreboard;
		for (int i = 0; i < min(7, static_cast<int>(algorithm.scoreboard.size())); ++i)
		{
			geneScoreboard += "#" + to_string(i + 1) + " : " + to_string(algorithm.scoreboard[i].score) + " {";
			for (unsigned int j = 0; j < algorithm.scoreboard[i].genes.size(); ++j)
			{
				geneScoreboard += to_string(algorithm.scoreboard[i].genes[j].data) + ", ";
			}
			geneScoreboard.resize(geneScoreboard.size() - 2);
			geneScoreboard += "}\n";
		}

		drawTextShadowed(fnt, win.width() - textOffset, textOffset, "FPS : " + to_string(framerateGet()), zoomRatio, zoomRatio, 2, 0);
		drawTextShadowed(fnt, textOffset, textOffset, "Generation #" + to_string(algorithm.getGeneration()) +
													  "\nTotal canonballs : " + to_string(algorithm.balls.size()) +
													  "\n\nLeading genes (last iteration) :\n" + geneScoreboard +
													  "\nCurrent score (last iteration) : " + to_string(algorithm.getLastIterationScore()) +
													  "\nScore to beat (last generation) : " + to_string(algorithm.getLastGenerationScore()), zoomRatio, zoomRatio);

		glColor4ub(0, 0, 0, 50);
		drawRectangle(win.width() / 2, win.height() - (zoomRatio * 128), win.width() - (16 * zoomRatio), win.height() - (zoomRatio * 16));
		glColor4ub(255, 0, 0, 255);

		for (float i = 0; i < std::max(1, static_cast<int>(algorithm.scoreevol.size())) - 1; ++i)
		{
			float dotpx =  (i   * (win.width() / 2.f - 32.f * zoomRatio)) / algorithm.scoreevol.size();
			float dotpx2 = (1.f * (win.width() / 2.f - 32.f * zoomRatio)) / algorithm.scoreevol.size();

			float x = win.width() / 2.f + (8.f * zoomRatio) + dotpx;

			float dotpy =  (max(19800, static_cast<int>(algorithm.scoreevol[i]))     - 19800.f) / 200.f * 96.f * zoomRatio;
			float dotpy2 = (max(19800, static_cast<int>(algorithm.scoreevol[i + 1])) - 19800.f) / 200.f * 96.f * zoomRatio;

			float ybase = win.height() - (zoomRatio * 24);

            drawLine(x, ybase - dotpy, x + dotpx2, ybase - dotpy2);
		}

		glColor4ub(255, 255, 255, 255);

        win.screenRefresh();
        framerateUpdate();
    }

    return 0;
}
