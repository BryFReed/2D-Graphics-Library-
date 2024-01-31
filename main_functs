/*
 *  Copyright 2023 <me>
 */

#include "include/GCanvas.h"
#include "include/GRect.h"
#include "include/GColor.h"
#include "include/GBitmap.h"
#include <list>
#include "include/GMatrix.h"
#include "include/GShader.h"
#include "include/GPath.h"
#include "include/GFinal.h"

std::unique_ptr<GShader> GCreateTrigShader(GPoint p1, GPoint p2, GPoint p3, GColor cols[], GPoint texs[]);
    
using namespace std;
GPixel bilerp(GPixel surr[], float x, float y){
    GPixel cent = surr[0];
    GPixel upp = surr[1];
    GPixel rightp = surr[2];
    GPixel down = surr[3];
    GPixel left = surr[4];
    GPixel ret;
    float x_diff;
    float y_diff;
    bool right;
    bool up;
    float r = 0;
    float g = 0;
    float b = 0;
    float a = 0;
    if(x <= 0.5){
        x_diff = (0.5 - x)/0.5;
        right = false;
    }
    else{
        x_diff = abs((0.5-x)/0.5);
        right = true;
    }
    if(y <= 0.5){
        y_diff = (0.5 - y)/0.5;
        up = false;
    }
    else{
        y_diff = abs((0.5-y)/0.5);
        up = true;
    }
    if(up){
        r += GPixel_GetR(upp) * y_diff;
        r += GPixel_GetR(cent) * (1-y_diff);
        g += GPixel_GetG(upp) * y_diff;
        g += GPixel_GetG(cent) * (1-y_diff);
        b += GPixel_GetB(upp) * y_diff;
        b += GPixel_GetB(cent) * (1-y_diff);
        a += GPixel_GetA(upp) * y_diff;
        a += GPixel_GetA(cent) * (1-y_diff);
    }
    else{
        r += GPixel_GetR(down) * y_diff;
        r += GPixel_GetR(cent) * (1-y_diff);
        g += GPixel_GetG(down) * y_diff;
        g += GPixel_GetG(cent) * (1-y_diff);
        b += GPixel_GetB(down) * y_diff;
        b += GPixel_GetB(cent) * (1-y_diff);
        a += GPixel_GetA(down) * y_diff;
        a += GPixel_GetA(cent) * (1-y_diff);
    }
    if(right){
        r += GPixel_GetR(rightp) * x_diff;
        r += GPixel_GetR(cent) * (1-x_diff);
        g += GPixel_GetG(rightp) * x_diff;
        g += GPixel_GetG(cent) * (1-x_diff);
        b += GPixel_GetB(rightp) * x_diff;
        b += GPixel_GetB(cent) * (1-x_diff);
        a += GPixel_GetA(rightp) * x_diff;
        a += GPixel_GetA(cent) * (1-x_diff);
    }
    else{
        r += GPixel_GetR(left) * x_diff;
        r += GPixel_GetR(cent) * (1-x_diff);
        g += GPixel_GetG(left) * x_diff;
        g += GPixel_GetG(cent) * (1-x_diff);
        b += GPixel_GetB(left) * x_diff;
        b += GPixel_GetB(cent) * (1-x_diff);
        a += GPixel_GetA(left) * x_diff;
        a += GPixel_GetA(cent) * (1-x_diff);
    }
    r *= 0.5;
    g *= 0.5;
    b *= 0.5;
    a *= 0.5;
    ret = GPixel_PackARGB(a, r,g,b);
    return ret;
}
class ProxyShader : public GShader {
    GShader* fRealShader;
    GMatrix  fExtraTransform;
public:
    ProxyShader(GShader* shader, const GMatrix extraTransform)
        : fRealShader(shader), fExtraTransform(extraTransform) {}

    bool isOpaque() override { return this->fRealShader->isOpaque(); }

    bool setContext(const GMatrix& ctm) override {
        return this->fRealShader->setContext(ctm * this->fExtraTransform);
    }
    
    void shadeRow(int x, int y, int count, GPixel row[]) override {
        this->fRealShader->shadeRow(x, y, count, row);
    }
};
class DubShader : public GShader {
    GShader* fRealShader;
    GShader* other;
    GMatrix  fExtraTransform;
public:
    DubShader(GShader* shader, GShader* other,  const GMatrix extraTransform)
        : fRealShader(shader), other(other), fExtraTransform(extraTransform) {}

    bool isOpaque() override { return this->fRealShader->isOpaque(); }

    bool setContext(const GMatrix& ctm) override {
        return(this->fRealShader->setContext(ctm) && this->other->setContext(ctm));
    }
    
    void shadeRow(int x, int y, int count, GPixel row[]) override {
        this->fRealShader->shadeRow(x, y, count, row);
        GPixel thing[count];
        this->other->shadeRow(x, y, count, thing);
        for(int i = 0; i < count; i++){
            row[i] = GPixel_PackARGB(GPixel_GetA(row[i]) * GPixel_GetA(thing[i])/255,
                                     GPixel_GetR(row[i]) * GPixel_GetR(thing[i])/255,
                                     GPixel_GetG(row[i]) * GPixel_GetG(thing[i])/255,
                                     GPixel_GetB(row[i]) * GPixel_GetB(thing[i])/255);
        }
    }
};

void GPath::addRect(const GRect& r, Direction dir){
    this->moveTo(r.left, r.top);
    switch(dir){
        case Direction::kCW_Direction:
            this->lineTo(r.right, r.top);
            this->lineTo(r.right, r.bottom);
            this->lineTo(r.left, r.bottom);
            break;
        case Direction::kCCW_Direction:
            this->lineTo(r.left, r.bottom);
            this->lineTo(r.right, r.bottom);
            this->lineTo(r.right, r.top);
            break;
    }
}
void GPath::addPolygon(const GPoint pts[], int count){
    this->moveTo(pts[0]);
    for(int i = 1; i < count; i++){
        this->lineTo(pts[i]);
    }
}
void GPath::addCircle(GPoint center, float radius, GPath::Direction dir) {
    if(dir == GPath::Direction::kCW_Direction){
        GPoint a = {1,0};
        GPoint b = {1, tan(M_PI_4 / 2)};
        GPoint c = {sqrt(2) / 2, sqrt(2) / 2};
        GMatrix mx = GMatrix::Translate(center.x, center.y) * GMatrix::Scale(radius, radius);
        this->moveTo(mx * a);
        for (int i = 0; i < 8; ++i) {
            this->quadTo(mx * (b), mx * (c));
            GMatrix rotate = mx.Rotate(M_PI_4);
            a = rotate * (a);
            b = rotate * (b);
            c = rotate * (c);
        }
    }
    else{
        GPoint a = {1,0};
        GPoint b = {1, -tan(M_PI_4 / 2)};
        GPoint c = {sqrt(2) / 2, -sqrt(2) / 2};
        GMatrix mx = GMatrix::Translate(center.x, center.y) * GMatrix::Scale(radius, radius);
        this->moveTo(mx * a);
        for (int i = 0; i < 8; ++i) {
            this->quadTo(mx * (b), mx * (c));
            GMatrix rotate = mx.Rotate(-M_PI_4);
            a = rotate * (a);
            b = rotate * (b);
            c = rotate * (c);
        }
        
    }
}
GRect GPath::bounds() const{
    if (fPts.empty()) {
                return {0.0, 0.0, 0.0, 0.0};
            }

            float minX = fPts[0].x;
            float minY = fPts[0].y;
            float maxX = fPts[0].x;
            float maxY = fPts[0].y;

            for (const auto& point : fPts) {
                minX = std::min(minX, point.x);
                minY = std::min(minY, point.y);
                maxX = std::max(maxX, point.x);
                maxY = std::max(maxY, point.y);
            }

            return {minX, minY, maxX, maxY};
}
void GPath::transform(const GMatrix& m){
    for(int i = 0; i < fPts.size(); i ++){
        fPts[i] = m*fPts[i];
    }
}
GPoint lerp(GPoint a, GPoint b, float t){
        return a + (b-a)*t;
}
void GPath::ChopQuadAt(const GPoint src[3], GPoint dst[5], float t){
    GPoint a = src[0];
    GPoint b = src[1];
    GPoint c = src[2];
    dst[0] = a;
    dst[1] = lerp(a, b, t);
    dst[3] = lerp(b, c, t);
    dst[4] = c;
    dst[2] = lerp(dst[1], dst[3], t);
}
void GPath::ChopCubicAt(const GPoint src[4], GPoint dst[7], float t){
    GPoint a = src[0];
    GPoint b = src[1];
    GPoint c = src[2];
    GPoint d = src[3];
    dst[0] = a;
    dst[1] = lerp(a, b, t);
    dst[2] = lerp(dst[1], lerp(b,c,t) ,t);
    dst[5] = lerp(c, d, t);
    dst[6] = d;
    dst[4] = lerp(lerp(b,c,t), dst[5], t);
    dst[3] = lerp(dst[2], dst[4], t);
}
struct GSeg {
    GSeg(GPoint topp, GPoint bott){
        double deltaY = topp.y - bott.y;
        double deltaX = topp.x - bott.x;
        m = deltaX/deltaY;
        b = topp.x - topp.y * m;
        this->top.x = GRoundToInt(topp.x);
        this->top.y = GRoundToInt(topp.y);

        this->bot.x = GRoundToInt(bott.x);
        this->bot.y = GRoundToInt(bott.y);

    }
    bool isValid(int y) const {
        
        return(y>= this->bot.y && y < this->top.y);
    }

    double m;
    double b;
    GPoint top;
    GPoint bot;
    int wind = 0;
};

GPath stroke(const GPoint points[], float width, int count, bool close){
    GPath nu;
    GPoint pts[2];
    vector<GPoint> lines;
    GPoint recs[4];
    GPoint flat[3];
    GPoint pointy[3];
    float yv;
    float xv;
    float mag;
    int other_count = count;
    if(close == false){
        other_count -=1;
        nu.addCircle(points[0], width/2, GPath::Direction::kCCW_Direction);
    }
    for(int i = 0; i < other_count; i++){
        pts[0] = points[i];
        pts[1] = points[(i+1)% count];
        yv = pts[1].x - pts[0].x;
        xv = pts[1].y - pts[0].y;
        mag = sqrt(xv * xv + yv * yv);
        xv /= mag; xv *= width/2;
        yv /= mag; yv *= width/2;
        recs[1] = {pts[0].x - xv, pts[0].y + yv};
        recs[0] = {pts[0].x + xv, pts[0].y - yv};
        recs[2] = {pts[1].x - xv, pts[1].y + yv};
        recs[3] = {pts[1].x + xv, pts[1].y - yv};
        nu.addPolygon(recs, 4);
        nu.addCircle(pts[1], width/2, GPath::Direction::kCCW_Direction);
    }
    return nu;
}

class RadialShader : public GShader {
    GPoint center;
    float radius;
    int count;
    GMatrix ctm;
    GShader::TileMode mode;
    vector<GColor> colors;

public:
    RadialShader(GPoint cent, float rad, const GColor cols[], int cnt, GShader::TileMode md)
        : center(cent), radius(rad), count(cnt), mode(md) {
            this->colors.resize(count+1);

            for(int i = 0; i < count; i++){
                this->colors[i] = cols[i];
            }
            this->colors[count] = cols[count-1];
        }

    bool isOpaque() override { return true; }

    bool setContext(const GMatrix& ctm) override {
        this->ctm = ctm;
        this->ctm.invert(&this->ctm);
        return true;
    }
    
    void shadeRow(int x, int y, int count, GPixel row[]) override {
        for(int i = 0; i < count; i++){
            GPoint p = ctm * GPoint{ x + 0.5, y + 0.5 };
            float xdis = center.x - p.x;
            float ydis = center.y - p.y;
            float dis = sqrt(xdis * xdis + ydis * ydis);
            bool skip = false;
            if(this->mode == GShader::TileMode::kClamp){
            }
          if(this->mode == GShader::TileMode::kRepeat){
                dis = fmod(dis, radius);
          }
            if(this->mode == GShader::TileMode::kMirror){
                int chec = dis/radius;
                if(chec%2 == 0){
                    dis = fmod(dis, radius);
                }
                else{
                    dis = radius - fmod(dis, radius);
                }
            }
    
            float temp = dis/radius;
            temp *= this->colors.size() -2;
            float t = temp - int(temp);
            GPixel packedPixel = GPixel_PackARGB(GRoundToInt(colors[int(temp)].a * (1 - t) * 255 + colors[(int(temp + 1))].a * t * 255)
                                                , GRoundToInt(colors[int(temp)].r * (1 - t) * colors[int(temp)].a * 255 + colors[(int(temp + 1))].r * colors[(int(temp + 1))].a * t * 255)
                                                , GRoundToInt(colors[int(temp)].g * (1 - t) * colors[int(temp)].a * 255 + colors[(int(temp + 1))].g * colors[(int(temp + 1))].a * t * 255)
                                                , GRoundToInt(colors[int(temp)].b * (1 - t) * colors[int(temp)].a * 255 + colors[(int(temp + 1))].b * colors[(int(temp + 1))].a * t * 255));

            row[i] = packedPixel;
            p.x += 1;
        }
    }
};
class myFinal : public GFinal{
public:
    myFinal(){
        
    }
    GPath strokePolygon(const GPoint pts[], int count, float width, bool isClosed){
        return stroke(pts, width, count, isClosed);
    }
    std::unique_ptr<GShader> createColorMatrixShader(const GColorMatrix&,
                                                     GShader* realShader) {
        return nullptr;
    }
    std::unique_ptr<GShader> createRadialGradient(GPoint center, float radius,
                                                  const GColor colors[], int count,
                                                  GShader::TileMode mode) {
        return unique_ptr<GShader>(new RadialShader(center, radius, colors, count, mode));
    }
    void drawQuadraticCoons(GCanvas*, const GPoint pts[8], const GPoint tex[4],
                            int level, const GPaint&) {}
};
std::unique_ptr<GFinal> GCreateFinal(){
    return std::unique_ptr<GFinal>(new myFinal());
}
bool compSeg(GSeg thing, GSeg other)
{
    if(thing.top.y == other.top.y){
        return(thing.bot.x < other.bot.x);
    }
    return (thing.top.y > other.top.y);
}

int findIntersection(const GSeg& segment, double y) {
    int x = segment.m*(y+.5) + segment.b;
    return x;
    }
bool compSeg1(GSeg thing, GSeg other)
{
    if(thing.bot.y == other.bot.y){
        return(findIntersection(thing, thing.bot.y) < findIntersection(other, other.bot.y));
    }
    return (thing.bot.y < other.bot.y);
}

static inline GPixel kClear(GPixel source, GPixel other){
    return 0;
}
static inline GPixel kSrc(GPixel source, GPixel other){
    return source;
}
static inline GPixel kDst(GPixel source, GPixel other){
    return other;
}
static inline GPixel kSrcOver(GPixel source, GPixel other){
    unsigned char newa, newr, newg, newb;
    newa = GPixel_GetA(source) + (255 - GPixel_GetA(source)) * GPixel_GetA(other) / 255;
    newr = GPixel_GetR(source) + (255 - GPixel_GetA(source)) * GPixel_GetR(other) / 255;
    newg = GPixel_GetG(source) + (255 - GPixel_GetA(source)) * GPixel_GetG(other) / 255;
    newb = GPixel_GetB(source) + (255 - GPixel_GetA(source)) * GPixel_GetB(other) / 255;
    GPixel packed = GPixel_PackARGB(newa , newr, newg, newb);
    return packed;
}
static inline GPixel kDstOver(GPixel source, GPixel other){
    unsigned char newa, newr, newg, newb;
    newa = GPixel_GetA(other) + (255 - GPixel_GetA(other)) * GPixel_GetA(source) / 255;
    newr = GPixel_GetR(other) + (255 - GPixel_GetA(other)) * GPixel_GetR(source) / 255;
    newg = GPixel_GetG(other) + (255 - GPixel_GetA(other)) * GPixel_GetG(source) / 255;
    newb = GPixel_GetB(other) + (255 - GPixel_GetA(other)) * GPixel_GetB(source) / 255;
    GPixel packed = GPixel_PackARGB(newa , newr, newg, newb);
    return packed;
}
static inline GPixel kSrcIn(GPixel source, GPixel other){
    unsigned char newa, newr, newg, newb;
    newa = GPixel_GetA(source) * GPixel_GetA(other) / 255;
    newr = GPixel_GetR(source) * GPixel_GetA(other) / 255;
    newg = GPixel_GetG(source) * GPixel_GetA(other) / 255;
    newb = GPixel_GetB(source) * GPixel_GetA(other) / 255;
    GPixel packed = GPixel_PackARGB(newa , newr, newg, newb);
    return packed;
}
static inline GPixel kDstIn(GPixel source, GPixel other){
    unsigned char newa, newr, newg, newb;
    newa = GPixel_GetA(other) * GPixel_GetA(source) / 255;
    newr = GPixel_GetR(other) * GPixel_GetA(source) / 255;
    newg = GPixel_GetG(other) * GPixel_GetA(source) / 255;
    newb = GPixel_GetB(other) * GPixel_GetA(source) / 255;
    GPixel packed = GPixel_PackARGB(newa , newr, newg, newb);
    return packed;
}
static inline GPixel kSrcOut(GPixel source, GPixel other){
    unsigned char newa, newr, newg, newb;
    newa = (255 - GPixel_GetA(other)) * GPixel_GetA(source) / 255;
    newr = (255 - GPixel_GetA(other)) * GPixel_GetR(source) / 255;
    newg = (255 - GPixel_GetA(other)) * GPixel_GetG(source) / 255;
    newb = (255 - GPixel_GetA(other)) * GPixel_GetB(source) / 255;
    GPixel packed = GPixel_PackARGB(newa , newr, newg, newb);
    return packed;
}
static inline GPixel kDstOut(GPixel source, GPixel other){
    unsigned char newa, newr, newg, newb;
    newa = (255 - GPixel_GetA(source)) * GPixel_GetA(other) / 255;
    newr = (255 - GPixel_GetA(source)) * GPixel_GetR(other) / 255;
    newg = (255 - GPixel_GetA(source)) * GPixel_GetG(other) / 255;
    newb = (255 - GPixel_GetA(source)) * GPixel_GetB(other) / 255;
    GPixel packed = GPixel_PackARGB(newa , newr, newg, newb);
    return packed;
}
static inline GPixel kSrcATop(GPixel source, GPixel other){
    unsigned char newa, newr, newg, newb;
    newa = GPixel_GetA(source) * GPixel_GetA(other) / 255 + (255 - GPixel_GetA(source)) * GPixel_GetA(other) / 255;
    newr = GPixel_GetR(source) * GPixel_GetA(other) / 255 + (255 - GPixel_GetA(source)) * GPixel_GetR(other) / 255;
    newg = GPixel_GetG(source) * GPixel_GetA(other) / 255 + (255 - GPixel_GetA(source)) * GPixel_GetG(other) / 255;
    newb = GPixel_GetB(source) * GPixel_GetA(other) / 255 + (255 - GPixel_GetA(source)) * GPixel_GetB(other) / 255;
    GPixel packed = GPixel_PackARGB(newa , newr, newg, newb);
    return packed;
}
static inline GPixel kDstATop(GPixel source, GPixel other){
    unsigned char newa, newr, newg, newb;
    newa = GPixel_GetA(other) * GPixel_GetA(source) / 255 + (255 - GPixel_GetA(other)) * GPixel_GetA(source) / 255;
    newr = GPixel_GetR(other) * GPixel_GetA(source) / 255 + (255 - GPixel_GetA(other)) * GPixel_GetR(source) / 255;
    newg = GPixel_GetG(other) * GPixel_GetA(source) / 255 + (255 - GPixel_GetA(other)) * GPixel_GetG(source) / 255;
    newb = GPixel_GetB(other) * GPixel_GetA(source) / 255 + (255 - GPixel_GetA(other)) * GPixel_GetB(source) / 255;
    GPixel packed = GPixel_PackARGB(newa , newr, newg, newb);
    return packed;
}
static inline GPixel kXor(GPixel source, GPixel other){
    unsigned char newa, newr, newg, newb;
    newa = (255 - GPixel_GetA(source)) * GPixel_GetA(other) / 255 + (255 - GPixel_GetA(other)) * GPixel_GetA(source) / 255;
    newr = (255 - GPixel_GetA(source)) * GPixel_GetR(other) / 255 + (255 - GPixel_GetA(other)) * GPixel_GetR(source) / 255;
    newg = (255 - GPixel_GetA(source)) * GPixel_GetG(other) / 255 + (255 - GPixel_GetA(other)) * GPixel_GetG(source) / 255;
    newb = (255 - GPixel_GetA(source)) * GPixel_GetB(other) / 255 + (255 - GPixel_GetA(other)) * GPixel_GetB(source) / 255;
    GPixel packed = GPixel_PackARGB(newa , newr, newg, newb);
    return packed;
}
template <typename BFunc> void blendFast(GPixel source, GPixel other[], int count, BFunc func){
    for(int i = 0; i < count; i++){
        other[i] = func(source, other[i]);
    }
}
static void blendArr(GPixel source, GPixel other[], int count, GBlendMode type){
    
    switch (type) {
        case GBlendMode::kClear:
            blendFast(source,other,count,kClear);
            break;

        case GBlendMode::kSrc:
            blendFast(source,other,count,kSrc);
            break;

        case GBlendMode::kDst:
            blendFast(source,other,count,kDst);
            break;

        case GBlendMode::kSrcOver:
            blendFast(source,other,count,kSrcOver);
            break;

        case GBlendMode::kDstOver:
            blendFast(source,other,count,kDstOver);
            break;

        case GBlendMode::kSrcIn:
            blendFast(source,other,count,kSrcIn);
            break;


        case GBlendMode::kDstIn:
            blendFast(source,other,count,kDstIn);
            break;


        case GBlendMode::kSrcOut:
            blendFast(source,other,count,kSrcOut);
            break;


        case GBlendMode::kDstOut:
            blendFast(source,other,count,kDstOut);
            break;


        case GBlendMode::kSrcATop:
            blendFast(source,other,count,kSrcATop);
            break;


        case GBlendMode::kDstATop:
            blendFast(source,other,count,kDstATop);
            break;


        case GBlendMode::kXor:
            blendFast(source,other,count,kXor);
            break;

    }
    
}


static GPixel blend(GPixel source, GPixel other, GBlendMode type){

        switch (type) {
            case GBlendMode::kClear:
                return kClear(source, other);

            case GBlendMode::kSrc:
                return kSrc(source, other);


            case GBlendMode::kDst:
                return kDst(source, other);

            case GBlendMode::kSrcOver:
                return kSrcOver(source, other);

            case GBlendMode::kDstOver:
                return kDstOver(source, other);


            case GBlendMode::kSrcIn:
                return kSrcIn(source, other);


            case GBlendMode::kDstIn:
                return kDstIn(source, other);


            case GBlendMode::kSrcOut:
                return kSrcOut(source, other);


            case GBlendMode::kDstOut:
                return kDstOut(source, other);


            case GBlendMode::kSrcATop:
                return kSrcATop(source, other);


            case GBlendMode::kDstATop:
                return kDstATop(source, other);


            case GBlendMode::kXor:
                return kXor(source, other);

        }
        
    
    
}
GBlendMode choseBlend(GPaint thing){
    float sa = thing.getColor().a;
    GBlendMode mode = thing.getBlendMode();
    GShader* s = thing.getShader();
    if(s){
        sa = 0.5;
        if (s->isOpaque()){
            sa = 1;
        }
    }
    if(sa == 0){
        if(mode == GBlendMode::kSrc || mode == GBlendMode::kDstIn || mode == GBlendMode::kSrcIn ){
            mode = GBlendMode::kClear;
        }
        if(mode == GBlendMode::kDstOut){
            mode = GBlendMode::kDst;
        }
        if(mode == GBlendMode::kDstATop){
            mode = GBlendMode::kClear;
        }
        if(mode == GBlendMode::kSrcOver){
            mode = GBlendMode::kDst;
        }
        if(mode == GBlendMode::kDstOver){
            mode = GBlendMode::kDst;
        }
        
    }
    if(sa == 1){
        if(mode == GBlendMode::kSrcOver){
            mode = GBlendMode::kSrc;
        }
        if(mode == GBlendMode::kDstIn){
            mode = GBlendMode::kDst;
        }
        if(mode == GBlendMode::kDstOut){
            mode = GBlendMode::kClear;
        }
        if(mode == GBlendMode::kSrcATop){
            mode = GBlendMode::kSrcIn;
        }
        if(mode == GBlendMode::kDstATop){
            mode = GBlendMode::kDstOver;
        }
        if(mode == GBlendMode::kXor){
            mode = GBlendMode::kSrcOut;
        }
    }
    return mode;
}
vector<GPoint> findQuad(GPoint a, GPoint b, GPoint c, int count){
    vector<GPoint> pts(count+1);
    int thing = 0;
    for(float i = 0; i < 1; i += (1.0f/count)){
        GPoint add = (1-i)*(a*(1-i) + 2*i*b)+ i*i*c;
        if(i == 0){
            add = a;
        }
        pts[thing] = add;
        thing += 1;
    }
    pts[count] = c;
    return pts;
}
vector<GPoint> findCubic(GPoint a, GPoint b, GPoint c, GPoint d, int count){
    vector<GPoint> pts(count+1);
    int thing = 0;
    for(float i = 0; i < 1; i+= (1.0f/count)){
        GPoint add = (1-i)*((1-i)*(1-i)*a+i*(3*(1-i)*b+3*i*c))+i*i*i*d;
        if(i == 0){
            add = a;
        }
        pts[thing] = add;
        thing += 1;
    }
    pts[count] = d;
    return pts;
}
class MyCanvas : public GCanvas {
    GMatrix ident = GMatrix(1, 0, 0, 0, 1, 0);
    list<GMatrix> stack;
public:
    MyCanvas(const GBitmap& device) : fDevice(device) {stack.push_back(ident);}
    
    void clear(const GColor& color) override {
        for(int i = 0; i < this->fDevice.width(); i++){
            for(int j = 0; j < this->fDevice.height(); j++){
                GPixel packedPixel = GPixel_PackARGB(color.a*255, color.r*255* color.a, color.g*255 * color.a, color.b*255 * color.a);
                fDevice.getAddr(i, j)[0] = packedPixel;
            }
                
        }
    }
    void drawPath(const GPath& path, const GPaint& color){
        if(color.getShader()){
            if(!color.getShader()->setContext(this->stack.back())){
                printf("imnot here/n");
                return;
            }
        }
        GBlendMode mode = choseBlend(color);
        if(mode == GBlendMode::kDst){
            return;
        }
        GPath::Edger e(path);
        GPoint pts[4];
        vector<GSeg> lines;
        GPath::Verb v;
        while((v = e.next(pts)) != GPath::Verb::kDone){
            if (v == GPath::Verb::kLine){
                bool add = true;
                pts[0] = stack.back() * pts[0];
                pts[1] = stack.back() * pts[1];
                if(pts[0].y > pts[1].y){
                    
                    GSeg thing(pts[0], pts[1]);
                    thing.wind = -1;
                    if(thing.top.y > this->fDevice.height()){
                        thing.top.y = this->fDevice.height();
                        if(thing.bot.y > this->fDevice.height()){
                            add = false;
                        }
                    }
                    
                    if(thing.bot.y < 0){
                        thing.bot.y = 0;
                        if(thing.top.y < 0){
                            add = false;
                        }
                    }
                    if(thing.top.x > this->fDevice.width()){
                        thing.top.x = this->fDevice.width();
                    }
                    if(thing.bot.x > this->fDevice.width()){
                        thing.bot.x = this->fDevice.width();
                    }
                    
                    if(thing.bot.y < 0){
                        thing.bot.y = 0;
                        if(thing.top.y < 0){
                            thing.top.y = 0;
                        }
                    }
                    if(thing.top.x < 0){
                        thing.top.x = 0;
                    }
                    if(thing.bot.x < 0){
                        thing.bot.x = 0;
                    }
                    if(add && thing.bot.y != thing.top.y){
                        lines.push_back(thing);
                    }
                }
                else{
                    GSeg thing(pts[1], pts[0]);
                    thing.wind = 1;
                    if(thing.top.y > this->fDevice.height()){
                        thing.top.y = this->fDevice.height();
                        if(thing.bot.y > this->fDevice.height()){
                            thing.bot.y = this->fDevice.height();
                        }
                    }
                    if(thing.top.x > this->fDevice.width()){
                        thing.top.x = this->fDevice.width();
                    }
                    if(thing.bot.x > this->fDevice.width()){
                        thing.bot.x = this->fDevice.width();
                    }
                    
                    if(thing.bot.y < 0){
                        thing.bot.y = 0;
                        if(thing.top.y < 0){
                            thing.top.y = 0;
                        }
                    }
                    if(thing.top.x < 0){
                        thing.top.x = 0;
                    }
                    if(thing.bot.x < 0){
                        thing.bot.x = 0;
                    }
                    
                    if(add && thing.bot.y != thing.top.y){
                        lines.push_back(thing);
                    }
                }
            }
            if(v == GPath::Verb::kQuad){
                stack.back().mapPoints(pts, 3);
            bool add = true;
            GPoint e = (pts[0] - 2*pts[1] + pts[2])* (0.25f);
            int num_segs = (int)ceil(sqrt((sqrt(e.x*e.x + e.y*e.y))/0.25f));
            vector<GPoint> qPoints = findQuad(pts[0], pts[1], pts[2], num_segs);
                for(int z = 0; z < qPoints.size() - 1; z++){
                GPoint one = qPoints[z];
                GPoint two = qPoints[z+1];
                if(one.y > two.y){
                    
                    GSeg thing(one, two);
                    thing.wind = -1;
                    if(thing.top.y > this->fDevice.height()){
                        thing.top.y = this->fDevice.height();
                        if(thing.bot.y > this->fDevice.height()){
                            add = false;
                        }
                    }
                    
                    if(thing.bot.y < 0){
                        thing.bot.y = 0;
                        if(thing.top.y < 0){
                            add = false;
                        }
                    }
                    if(thing.top.x > this->fDevice.width()){
                        thing.top.x = this->fDevice.width();
                    }
                    if(thing.bot.x > this->fDevice.width()){
                        thing.bot.x = this->fDevice.width();
                    }
                    
                    if(thing.bot.y < 0){
                        thing.bot.y = 0;
                        if(thing.top.y < 0){
                            thing.top.y = 0;
                        }
                    }
                    if(thing.top.x < 0){
                        thing.top.x = 0;
                    }
                    if(thing.bot.x < 0){
                        thing.bot.x = 0;
                    }
                    if(add && thing.bot.y != thing.top.y){
                        lines.push_back(thing);
                    }
                }
                else{
                    GSeg thing(two, one);
                    thing.wind = 1;
                    if(thing.top.y > this->fDevice.height()){
                        thing.top.y = this->fDevice.height();
                        if(thing.bot.y > this->fDevice.height()){
                            thing.bot.y = this->fDevice.height();
                        }
                    }
                    if(thing.top.x > this->fDevice.width()){
                        thing.top.x = this->fDevice.width();
                    }
                    if(thing.bot.x > this->fDevice.width()){
                        thing.bot.x = this->fDevice.width();
                    }
                    
                    if(thing.bot.y < 0){
                        thing.bot.y = 0;
                        if(thing.top.y < 0){
                            thing.top.y = 0;
                        }
                    }
                    if(thing.top.x < 0){
                        thing.top.x = 0;
                    }
                    if(thing.bot.x < 0){
                        thing.bot.x = 0;
                    }
                    
                    if(add && thing.bot.y != thing.top.y){
                        lines.push_back(thing);
                    }
                }
            }
        }
            if(v == GPath::Verb::kCubic){
                stack.back().mapPoints(pts,4);
                bool add = true;
                GPoint E;
                GPoint E0 = pts[0] - 2*pts[1] + pts[2];
                GPoint E1 = pts[1] - 2*pts[2] + pts[3];
                E.x = max(abs(E0.x), abs(E1.x));
                E.y = max(abs(E0.y), abs(E1.y));
                int num_segs = (int)ceil(sqrt((3*sqrt(E.x*E.x + E.y*E.y))));
                vector<GPoint> qPoints = findCubic(pts[0], pts[1], pts[2], pts[3], num_segs);
                for(int z = 0; z < qPoints.size() - 1; z++){
                    GPoint one =qPoints[z];
                    GPoint two = qPoints[z+1];
                    if(one.y > two.y){
                        
                        GSeg thing(one, two);
                        thing.wind = -1;
                        if(thing.top.y > this->fDevice.height()){
                            thing.top.y = this->fDevice.height();
                            if(thing.bot.y > this->fDevice.height()){
                                add = false;
                            }
                        }
                        
                        if(thing.bot.y < 0){
                            thing.bot.y = 0;
                            if(thing.top.y < 0){
                                add = false;
                            }
                        }
                        if(thing.top.x > this->fDevice.width()){
                            thing.top.x = this->fDevice.width();
                        }
                        if(thing.bot.x > this->fDevice.width()){
                            thing.bot.x = this->fDevice.width();
                        }
                        
                        if(thing.bot.y < 0){
                            thing.bot.y = 0;
                            if(thing.top.y < 0){
                                thing.top.y = 0;
                            }
                        }
                        if(thing.top.x < 0){
                            thing.top.x = 0;
                        }
                        if(thing.bot.x < 0){
                            thing.bot.x = 0;
                        }
                        if(add && thing.bot.y != thing.top.y){
                            lines.push_back(thing);
                        }
                    }
                    else{
                        GSeg thing(two, one);
                        thing.wind = 1;
                        if(thing.top.y > this->fDevice.height()){
                            thing.top.y = this->fDevice.height();
                            if(thing.bot.y > this->fDevice.height()){
                                thing.bot.y = this->fDevice.height();
                            }
                        }
                        if(thing.top.x > this->fDevice.width()){
                            thing.top.x = this->fDevice.width();
                        }
                        if(thing.bot.x > this->fDevice.width()){
                            thing.bot.x = this->fDevice.width();
                        }
                        
                        if(thing.bot.y < 0){
                            thing.bot.y = 0;
                            if(thing.top.y < 0){
                                thing.top.y = 0;
                            }
                        }
                        if(thing.top.x < 0){
                            thing.top.x = 0;
                        }
                        if(thing.bot.x < 0){
                            thing.bot.x = 0;
                        }
                        
                        if(add && thing.bot.y != thing.top.y){
                            lines.push_back(thing);
                        }
                    }
                }
            }
            
        }
        if(lines.size() < 2){
            return;
        }
        int ymax =0;
        for(auto& L: lines){
            ymax = max(ymax, int(L.top.y+.5));
        }
        size_t count = lines.size();
        sort(lines.begin(), lines.begin() + lines.size(), compSeg1);
        for (int y = lines[0].bot.y; y < ymax; y++) {
            size_t i = 0;
            int w = 0;
            int L;
            while (i < lines.size() && lines[i].isValid(y)) {
                    int x = int(findIntersection(lines[i], y));
                if(x > this->fDevice.width()){
                    x = this->fDevice.width();
                }
                if(x < 0){
                    x = 0;
                }

                    if (w == 0) {
                        L = x;
                    }
                    w += lines[i].wind;  // +1 or -1
                    if (w == 0 && x > L) {
                        
                        GPixel packedPixel;
                        if(color.getShader()){
                            for(;L < x; L++){
                                color.getShader()->shadeRow(L, y, 1, &packedPixel);
                                fDevice.getAddr(L, y)[0] = blend(packedPixel, fDevice.getAddr(L, y)[0], mode);
                            }
                        }
                        else{
                            packedPixel = GPixel_PackARGB(color.getColor().a * 255, color.getColor().r * color.getColor().a * 255, color.getColor().g* color.getColor().a * 255, color.getColor().b * color.getColor().a * 255);
                            
                            blendArr(packedPixel,fDevice.getAddr(L, y), x-L, mode );
                        }
                    }

                    if (lines[i].isValid(y+1)) {
                        // if you track “currX”, bump it by slope here
                        i += 1;
                    } else {
                        assert(i < lines.size() && i >= 0);
                        lines.erase(lines.begin()+ i);
                        
                        // we’re done with this edge
                    }
            }
            while (i < lines.size() && lines[i].isValid(y+1)) {
                    i += 1;
                }

            sort(lines.begin(), lines.begin() + i, [y](GSeg thing, GSeg other){
                return (findIntersection(thing, y+1) < findIntersection(other, y+1));
            });
        }
    }
    

    void drawRect(const GRect& rect, const GPaint& color) override {
        GPoint one = GPoint{rect.left, rect.top};
        GPoint two = GPoint{rect.left, rect.bottom};
        GPoint three = GPoint{rect.right, rect.bottom};
        GPoint four = GPoint{rect.right, rect.top};
        GPoint ps[] = {one, two, three, four};
        if(color.getShader()){
            drawConvexPolygon(ps, 4, color);
        }
        else{
            
            
            GBlendMode mode = choseBlend(color);
            if(mode == GBlendMode::kDst){
                return;
            }
            int top = GRoundToInt(rect.y());
            if(top < 0){
                top = 0;
            }
            int left = GRoundToInt(rect.x());
            if(left < 0){
                left = 0;
            }
            int right = GRoundToInt(rect.right);
            if(right > this->fDevice.width()){
                right = this->fDevice.width();
            }
            
            int bot = GRoundToInt(rect.bottom);
            if(bot > this->fDevice.height()){
                bot = this->fDevice.height();
            }
            
            
            if(true){
                for(int i = top; i < bot; i++){
                    if(right > left){
                        GPixel packedPixel;
                        if(color.getShader()){
                            for(;left < right; left++){
                                color.getShader()->shadeRow(left, i, 1, &packedPixel);
                                fDevice.getAddr(left, i)[0] = blend(packedPixel, fDevice.getAddr(left, i)[0], mode);
                            }
                        }
                        else{
                            packedPixel = GPixel_PackARGB(color.getColor().a * 255, color.getColor().r * color.getColor().a * 255, color.getColor().g* color.getColor().a * 255, color.getColor().b * color.getColor().a * 255);
                            
                            blendArr(packedPixel,fDevice.getAddr(left, i), right-left, mode );
                        }
                    }
                }
            }
            else{
                for(int i = top; i < bot; i++){
                    for(int j = left; j < right; j++){
                        GPixel packedPixel = GPixel_PackARGB(color.getColor().a * 255, color.getColor().r * color.getColor().a * 255, color.getColor().g* color.getColor().a * 255, color.getColor().b * color.getColor().a * 255);
                        fDevice.getAddr(j, i)[0] = packedPixel;
                    }
                }
                
            }
        }

    }
    
    void drawConvexPolygon(const GPoint arry[], int count, const GPaint& color){
        if(color.getShader()){
                    if(!color.getShader()->setContext(this->stack.back())){
                        color.getShader()->setContext(this->stack.back());
                        printf("imnot here/n");
                        return;
                    }
                }
                if(count < 3){
                    return;
                }
                GBlendMode mode = choseBlend(color);
            if(mode == GBlendMode::kDst){
                return;
            }
                vector<GPoint> array;
                for(int i = 0; i < count; i++){
                    GPoint nu = stack.back()*arry[i];
                    array.push_back(nu);
                }
                vector<GSeg> lines;
                for(int i = 0; i < count; i++){
                    bool add = true;
                    if(array[i].y > array[(i+1) % count].y){
                        GSeg thing(array[i], array[(i+1) % count]);
                        if(thing.top.y > this->fDevice.height()){
                            thing.top.y = this->fDevice.height();
                            if(thing.bot.y > this->fDevice.height()){
                                add = false;
                            }
                        }

                        if(thing.bot.y < 0){
                            thing.bot.y = 0;
                            if(thing.top.y < 0){
                                add = false;
                            }
                        }
                        if(add && thing.bot.y != thing.top.y){
                            lines.push_back(thing);
                        }
                    }
                    else{
                        GSeg thing(array[(i+1) % count], array[i]);
                        if(thing.top.y > this->fDevice.height()){
                            thing.top.y = this->fDevice.height();
                            if(thing.bot.y > this->fDevice.height()){
                                thing.bot.y = this->fDevice.height();
                            }
                        }
                        if(thing.top.x > this->fDevice.width()){
                            thing.top.x = this->fDevice.width();
                        }
                        if(thing.bot.x > this->fDevice.width()){
                            thing.bot.x = this->fDevice.width();
                        }
                        
                        if(thing.bot.y < 0){
                            thing.bot.y = 0;
                            if(thing.top.y < 0){
                                thing.top.y = 0;
                            }
                        }
                        if(thing.top.x < 0){
                            thing.top.x = 0;
                        }
                        if(thing.bot.x < 0){
                            thing.bot.x = 0;
                        }
                        
                        if(add && thing.bot.y != thing.top.y){
                            lines.push_back(thing);
                        }
                    }
                }
                if(lines.size() < 2){
                    return;
                }
                sort(lines.begin(), lines.begin() + lines.size(), compSeg);
                GSeg left = lines[1];
                GSeg right = lines[0];
                if(lines[0].bot.x < lines[1].bot.x){
                     left = lines[0];
                     right = lines[1];
                }

                int lt = 2;
                
                for(int y = lines.front().top.y -1 ; y >= lines.back().bot.y; y--){
                    if(left.bot.y > y){
                        left = lines[lt];
                        lt += 1;
                    }
                    if(right.bot.y > y){
                        right = lines[lt];
                        lt += 1;
                    }
                    int x = findIntersection(left, y) + .5;
                    int x_r = findIntersection(right, y) + .5;
                    if(x < 0){
                        x = 0;
                    }
                    if(x_r < 0){
                        x_r = 0;
                    }
                    if(x_r > this->fDevice.width()){
                        x_r = this->fDevice.width();
                    }
                    if(x > this->fDevice.width()){
                        x = this->fDevice.width();
                    }
                    if(x > x_r){
                        int temp = x;
                        x = x_r;
                        x_r = temp;
                    }
                    assert(x <= x_r);
                    GPixel packedPixel;
                    if(x_r > x){
                        if(color.getShader()){
                            for(;x < x_r; x++){
                                color.getShader()->shadeRow(x, y, 1, &packedPixel);
                                fDevice.getAddr(x, y)[0] = blend(packedPixel, fDevice.getAddr(x, y)[0], mode);
                            }
                        }
                        else{
                            packedPixel = GPixel_PackARGB(color.getColor().a * 255, color.getColor().r * color.getColor().a * 255, color.getColor().g* color.getColor().a * 255, color.getColor().b * color.getColor().a * 255);
                            
                            blendArr(packedPixel,fDevice.getAddr(x, y), x_r-x, mode );
                        }
                    }
                       
                }
        
    }
    void drawMesh(const GPoint verts[], const GColor colors[], const GPoint texs[],int count, const int indices[], const GPaint& color){
        GPoint p1;
        GPoint p2;
        GPoint p3;
        GColor c1;
        GColor c2;
        GColor c3;
        GPoint t1;
        GPoint t2;
        GPoint t3;
        int n = 0;
        for(int i = 0; i < count; i++){
            GPoint pts[3];
            p1 = verts[indices[n]];
            p2 = verts[indices[n+1]];
            p3 = verts[indices[n+2]];
            pts[0] = p1; pts[1] = p2; pts[2] =p3;
            GColor cols[3];
            if(colors){
                c1 = colors[indices[n]];
                c2 = colors[indices[n+1]];
                c3 = colors[indices[n+2]];
                cols[0] = c1; cols[1] = c2; cols[2] = c3;
            }
            GPoint tex[3];
            if(texs){
                t1 = texs[indices[n]];
                t2 = texs[indices[n+1]];
                t3 = texs[indices[n+2]];
                tex[0] = t1; tex[1] = t2; tex[2] = t3;
            }
            
            GMatrix t  = GMatrix(t2.x - t1.x, t3.x - t1.x, t1.x,
                                 t2.y - t1.y, t3.y - t1.y, t1.y);
            GMatrix pp = GMatrix(p2.x - p1.x, p3.x - p1.x, p1.x,
                                 p2.y - p1.y, p3.y - p1.y, p1.y);
            t.invert(&t);
            auto trig = GCreateTrigShader(p1, p2, p3, cols, tex);
            auto proxy = ProxyShader(color.getShader(), pp * t);
            auto dub = DubShader(trig.get(), &proxy, pp * t);
            GPaint paint;
            if(colors && texs == nullptr){
                paint.setShader(trig.get());
            }
            else if(texs && colors == nullptr){
                paint.setShader(&proxy);

            }
            else{
                paint.setShader(&dub);
            }
            drawConvexPolygon(pts, 3, paint);
            n += 3;
            
        }
    }
    template <typename T> T quadlerp(float u, float v, T a, T b, T c, T d){
        return (1-v) * ((1-u)* a + u * b) + v*((1-u)*d + u*c);
    }
    void drawQuad(const GPoint verts[4], const GColor colors[4], const GPoint texs[4], int level, const GPaint& color){
        float u;
        float v;
        float space = 1.0f/(1 + level);
        for(int i = 0; i <= level; i++){
            for(int j = 0; j <= level; j++){
                u = i * space;
                v = j * space;
                GPoint pts[4];
                GPoint p1 = quadlerp(u, v, verts[0], verts[1], verts[2], verts[3]); pts[0] = p1;
                GPoint p2 = quadlerp(u+space, v, verts[0], verts[1], verts[2], verts[3]); pts[1] = p2;
                GPoint p3 = quadlerp(u+space, v+space, verts[0], verts[1], verts[2], verts[3]); pts[2] = p3;
                GPoint p4 = quadlerp(u, v+space, verts[0], verts[1], verts[2], verts[3]); pts[3] = p4;
                GColor colt[4];
                GPoint text[4];
                if(colors){
                    GColor c1 = quadlerp(u, v, colors[0], colors[1], colors[2], colors[3]); colt[0] = c1;
                    GColor c2 = quadlerp(u+space, v, colors[0], colors[1], colors[2], colors[3]); colt[1] = c2;
                    GColor c3 = quadlerp(u+space, v+space, colors[0], colors[1], colors[2], colors[3]); colt[2] = c3;
                    GColor c4 = quadlerp(u, v+space, colors[0], colors[1], colors[2], colors[3]); colt[3] = c4;
                    
                }
                if(texs){
                    GPoint t1 = quadlerp(u, v, texs[0], texs[1], texs[2], texs[3]); text[0] = t1;
                    GPoint t2 = quadlerp(u+space, v, texs[0], texs[1], texs[2], texs[3]); text[1] = t2;
                    GPoint t3 = quadlerp(u+space, v+space, texs[0], texs[1], texs[2], texs[3]); text[2] = t3;
                    GPoint t4 = quadlerp(u, v+space, texs[0], texs[1], texs[2], texs[3]); text[3] = t4;
                }
                int indc[] = {0,1,3};
                drawMesh(pts, colors ? colt : nullptr, texs ? text : nullptr, 1, indc, color);
                indc[0] = 2;
                drawMesh(pts, colors ? colt : nullptr, texs ? text : nullptr, 1, indc, color);

            }
        }
    }
    void save(){
        stack.push_back(stack.back());

    }
    void restore(){
        stack.pop_back();
    }
    void concat(const GMatrix& matrix){
        stack.back() = stack.back() * matrix;

    }
               
    

private:
    // Note: we store a copy of the bitmap
    const GBitmap fDevice;
        
};

std::unique_ptr<GCanvas> GCreateCanvas(const GBitmap& device) {
    return std::unique_ptr<GCanvas>(new MyCanvas(device));
}

std::string GDrawSomething(GCanvas* canvas, GISize dim) {
    GPath path;
    path.moveTo(0,0);
    path.lineTo(10,0);
    path.lineTo(0,10);
    canvas->drawPath(path, GPaint());
    

    return "tears in rain";
}

    GMatrix GMatrix::Translate(float tx, float ty) {
        return GMatrix(1, 0, tx, 0, 1, ty);
    }
GMatrix GMatrix::Scale(float sx, float sy) {
        return GMatrix(sx, 0, 0, 0, sy, 0);
    }
 GMatrix GMatrix::Rotate(float radians) {
        float cosVal = std::cos(radians);
        float sinVal = std::sin(radians);
        return GMatrix(cosVal, -sinVal, 0, sinVal, cosVal, 0);
    }
    
    /**
     *  Return the product of two matrices: a * b
     */
 GMatrix GMatrix::Concat(const GMatrix& a, const GMatrix& b){
        return GMatrix(
                       a[0] * b[0] + a[1] * b[3], a[0] * b[1] + a[1] * b[4], a[0] * b[2] + a[1] * b[5] + a[2],
                       a[3] * b[0] + a[4] * b[3], a[3] * b[1] + a[4] * b[4], a[3] * b[2] + a[4] * b[5] + a[5]
                       );
    }
    
    /*
     *  Compute the inverse of this matrix, and store it in the "inverse" parameter, being
     *  careful to handle the case where 'inverse' might alias this matrix.
     *
     *  If this matrix is invertible, return true. If not, return false, and ignore the
     *  'inverse' parameter.
     */
bool GMatrix::invert(GMatrix* inverse) const{
    float det = fMat[0] * (fMat[4] * 1 - fMat[5] * 0)
    - fMat[1] * (fMat[3] * 1 - fMat[5] * 0)
    + fMat[2] * (fMat[3] * 0 - fMat[4] * 0);
        if (det == 0) {
            return false;  // Matrix is not invertible
        }

        float invDet = 1.0f / det;
    float zero = fMat[0];
    inverse->fMat[0] = fMat[4] * 1 - fMat[5] * 0;
    float one = fMat[1];
    inverse->fMat[1] = fMat[2] * 0 - fMat[1] * 1;
    float two = fMat[2];
    inverse->fMat[2] = one * fMat[5] - fMat[2] * fMat[4];
    float three = fMat[3];
    inverse->fMat[3] = fMat[5] * 0 - three * 1;
    inverse->fMat[4] = zero * 1 - two * 0;
    inverse->fMat[5] = two * three - zero * fMat[5];
    for (int i = 0; i < 6; ++i) {
                inverse->fMat[i] = inverse->fMat[i] * invDet;
            }

    return true;
    }
    
    /**
     *  Transform the set of points in src, storing the resulting points in dst, by applying this
     *  matrix. It is the caller's responsibility to allocate dst to be at least as large as src.
     *
     *  [ a  b  c ] [ x ]     x' = ax + by + c
     *  [ d  e  f ] [ y ]     y' = dx + ey + f
     *  [ 0  0  1 ] [ 1 ]
     *
     *  Note: It is legal for src and dst to point to the same memory (however, they may not
     *  partially overlap). Thus the following is supported.
     *
     *  GPoint pts[] = { ... };
     *  matrix.mapPoints(pts, pts, count);
     */
void GMatrix::mapPoints(GPoint dst[], const GPoint src[], int count) const{
        for (int i = 0; i < count; ++i) {
            float x = src[i].x;
            float y = src[i].y;
            dst[i].x = fMat[0] * x + fMat[1] * y + fMat[2];
            dst[i].y = fMat[3] * x + fMat[4] * y + fMat[5];
        }
    }
GMatrix::GMatrix(){
    fMat[0] = 1;    fMat[1] = 0;    fMat[2] = 0;
    fMat[3] = 0;    fMat[4] = 1;    fMat[5] = 0;
}

class TrigShader: public GShader{
    GMatrix ctm;
    GMatrix m;
    GMatrix inv;
    GPoint u;
    GPoint v;
    vector<GColor> cols;
    vector<GPoint> texs;
public:
    TrigShader(GPoint p1, GPoint p2, GPoint p3, GColor cols[], GPoint texs[]){
        this->u = p2 - p1;
        this->v = p3 - p1;
        this->cols.resize(3);
        for(int i = 0; i < 3; i++){
            this->cols[i] = cols[i];
        }
        this->texs.resize(3);
        for(int i = 0; i < 3; i++){
            this->texs[i] = texs[i];
        }
        this->m = GMatrix(u.x, v.x, p1.x, u.y, v.y, p1.y);
    }
    bool setContext(const GMatrix& ctm){
        this->ctm = ctm;
        (this->ctm * this->m).invert(&this->inv);
        return true;
    }
    bool isOpaque(){
        return true;
    }
    void shadeRow(int x, int y, int count, GPixel out[]){
        for(int i = 0; i < count; i++){
            GPoint p = {x + 0.5, y + 0.5};
            p = this->inv * p;
            GColor c = p.x * cols[1] + p.y * cols[2] + (1 - p.x - p.y) * cols[0];
            if(c.a > 1){
                c.a = 1;
            }
            if(c.b > 1){
                c.b = 1;
            }
            if(c.g > 1){
                c.g = 1;
            }
            if(c.r > 1){
                c.r = 1;
            }
            GPixel pack = GPixel_PackARGB(c.a * 255, c.r * c.a * 255, c.g* c.a * 255, c.b * c.a * 255);
            out[i] = pack;
            x += 1;
        }
    }
};
std::unique_ptr<GShader> GCreateTrigShader(GPoint p1, GPoint p2, GPoint p3, GColor cols[], GPoint texs[]){
    return std::unique_ptr<GShader>(new TrigShader(p1, p2, p3, cols, texs));
}
class BitmapShader: public GShader{
    
    GBitmap map;
    GMatrix loc;
    GMatrix ctm;
    GShader::TileMode mode;
    float heightDiv;
    float widthDiv;
public:
    BitmapShader(GBitmap map, GMatrix loc, GShader::TileMode mode){
        this->map = map;
        this->loc = loc;
        this->mode = mode;
        this->heightDiv = 1.0f/(map.height());
        this->widthDiv = 1.0f/(map.width());
    }
    bool setContext(const GMatrix& ctm){
        GMatrix inv;
        bool thing = ctm.invert(&inv);
        this->ctm = inv;
        return thing;
    }
    virtual bool isOpaque(){
        return map.isOpaque();
    }
    void shadeRow(int x, int y, int count, GPixel out[]) {
        GPoint p = loc * ctm * GPoint{ x + 0.5, y + 0.5 };
        if(this->mode == GShader::TileMode::kClamp){
            for (int i = 0; i < count; ++i) {
                float ix = p.x;
                float iy = p.y;
                if(ix < 0){
                    ix = 0;
                }
                if(iy < 0){
                    iy = 0;
                }
                if(ix >= map.width()){
                    ix = map.width() -1;
                }
                if(iy >= map.height()){
                    iy = map.height() -1;
                }
                out[i] = *map.getAddr(ix, iy);
                p.x += (loc * ctm)[0]; // this is A
                p.y += (loc * ctm)[3]; // this is D
            }
        }
        if(this->mode == GShader::TileMode::kRepeat){
            for (int i = 0; i < count; ++i) {
                float ix = p.x;
                float iy = p.y;
                if(ix < 0){
                    ix *= this->widthDiv;
                    ix -= floor(ix);
                    ix = ix * (map.width()) + (map.width());
                }
                if(iy < 0){
                    iy *= this->heightDiv;
                    iy -= floor(iy);
                    iy = iy * (map.height()) + (map.height());
                }
                if(ix >= map.width()){
                    ix *= this->widthDiv;
                    ix -= floor(ix);
                    ix = ix * (map.width());
                }
                if(iy >= map.height()){
                    iy *= this->heightDiv;
                    iy -= floor(iy);
                    iy = iy * (map.height());
                }
                if(ix >= map.width()){
                    ix = map.width() -1;
                }
                if(iy >= map.height()){
                    iy = map.height() -1;
                }
                out[i] = *map.getAddr(ix, iy);
                p.x += (loc * ctm)[0]; // this is A
                p.y += (loc * ctm)[3]; // this is D
            }
        }
        if(this->mode == GShader::TileMode::kMirror){
            for (int i = 0; i < count; ++i) {
                float ix = p.x;
                float iy = p.y;
                if(ix < 0){
                    ix *= this->widthDiv;
                    if(int(ix) & 1){
                        ix -= floor(ix);
                        ix = abs(ix * (map.width()));
                    }
                    else{
                        ix -= floor(ix);
                        ix = ix * (map.width()) + (map.width());
                    }
                }
                if(iy < 0){
                    iy *= this->heightDiv;
                    if(int(iy) & 1){
                        iy -= floor(iy);
                        iy = abs(iy * (map.height()));
                    }
                    else{
                        iy -= floor(iy);
                        iy = iy * (map.height()) + (map.height());
                    }
                }
                if(ix >= map.width()){
                    ix *= this->widthDiv;
                    if(int(ix) & 1){
                        ix -= floor(ix);
                        ix = (map.width()) - ix * (map.width());
                    }
                    else{
                        ix -= floor(ix);
                        ix = ix * (map.width()-1);
                    }
                }
                if(iy >= map.height()){
                    iy *= this->heightDiv;
                    if(int(iy) & 1){
                        iy -= floor(iy);
                        iy = (map.height()) - iy * (map.height());
                    }
                    else{
                        iy -= floor(iy);
                        iy = iy * (map.height());
                    }
                }
                if(ix >= map.width()){
                    ix = map.width() -1;
                }
                if(iy >= map.height()){
                    iy = map.height() -1;
                }
                out[i] = *map.getAddr(ix, iy);
                p.x += (loc * ctm)[0]; // this is A
                p.y += (loc * ctm)[3]; // this is D
            }
        }
    }
};
std::unique_ptr<GShader> GCreateBitmapShader(const GBitmap& map, const GMatrix& localInverse, GShader::TileMode mode){
    return std::unique_ptr<GShader>(new BitmapShader(map, localInverse, mode));
}




class LinearShader : public GShader {

    GMatrix ctm;
    GMatrix ctmI;
    GMatrix thingI;
    GMatrix thing;
    vector<GColor> cols;
    GShader::TileMode mode;
    float widthDiv;
    
    
public:
    LinearShader(GPoint p0, GPoint p1, const GColor cols[], int count, GShader::TileMode mode){
        float dx = p1.x -p0.x;
        float dy = p1.y - p0.y;
        this->thing = GMatrix(dx, -dy, p0.x,
                              dy, dx, p0.y);
        this->thing.invert(&this->thingI);
        this->cols.resize(count+1);
        for(int i = 0; i < count; i++){
            this->cols[i] = cols[i];
        }
        this->cols[count] = this->cols[count-1];
        this->mode = mode;
        this->widthDiv = 1;
    }
        
        
        bool isOpaque() {
            return false;
        }
        
        bool setContext(const GMatrix& ctm)  {
            this->ctm = ctm;
            this->ctm.invert(&this->ctmI);
            return true;
        }
        
        void shadeRow(int x, int y, int count, GPixel out[])  {
            GPoint p = this->thingI* ctmI * GPoint{ x + 0.5, y + 0.5 } * (this->cols.size() -2);
            GMatrix gu = this->thingI * ctmI;
            for (int i = 0; i < count; ++i) {
                p = gu * GPoint{ x +i+ 0.5, y + 0.5 } ;
                float temp = p.x;
                if(temp < 0){
                    if(mode == GShader::TileMode::kClamp){
                        temp = 0;
                    }
                    if(mode == GShader::TileMode::kRepeat){
                        temp -= floor(temp);
                        temp = temp + 1;
                    }
                    if(mode == GShader::TileMode::kMirror){
                        temp -= floor(temp);
                        temp = abs(temp);
                    }
                }
                if(temp > (1)){
                    if(mode == GShader::TileMode::kClamp){
                        temp = 1;
                    }
                    if(mode == GShader::TileMode::kRepeat){
                        temp -= floor(temp);
                    }
                    if(mode == GShader::TileMode::kMirror){
                        temp -= floor(temp);
                        temp = 1-temp;
                    }
                }
                temp *= this->cols.size() -2;
                float t = temp - int(temp);
                GPixel packedPixel = GPixel_PackARGB(GRoundToInt(cols[int(temp)].a * (1 - t) * 255 + cols[(int(temp + 1))].a * t * 255)
                                                    , GRoundToInt(cols[int(temp)].r * (1 - t) * cols[int(temp)].a * 255 + cols[(int(temp + 1))].r * cols[(int(temp + 1))].a * t * 255)
                                                    , GRoundToInt(cols[int(temp)].g * (1 - t) * cols[int(temp)].a * 255 + cols[(int(temp + 1))].g * cols[(int(temp + 1))].a * t * 255)
                                                    , GRoundToInt(cols[int(temp)].b * (1 - t) * cols[int(temp)].a * 255 + cols[(int(temp + 1))].b * cols[(int(temp + 1))].a * t * 255));

                out[i] = packedPixel;
                p.x += (gu)[0]* (this->cols.size() -2); // this is A
            }
        }
        
    };

    
std::unique_ptr<GShader> GCreateLinearGradient(GPoint p0, GPoint p1, const GColor col[] , int count, GShader::TileMode mode)
 {
    return std::unique_ptr<GShader>(new LinearShader(p0,p1,col,count, mode ));
    }
